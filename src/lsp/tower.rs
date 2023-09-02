use dashmap::DashMap;
use lammps_analyser::error_finder::ErrorFinder;
use lammps_analyser::identifinder::IdentiFinder;
///! Alternative implementation for the lsp using the `tower-lsp` crate.
use tower_lsp::jsonrpc::Result;
use tower_lsp::lsp_types::*;
use tower_lsp::{Client, LanguageServer, LspService, Server};
use tree_sitter::{Parser, Tree};

/// Core LSP Application
/// TODO:
/// - [ ] Add Symbols
/// - [ ] Add Scemantic Tokens
///     - [ ] Create a token map using identifinder
///     - [ ] Add convinence methods for converting from this map to an LSP Scemantic Token
/// - [ ] Add Goto Definitions
#[derive(Debug)]
struct Backend {
    client: Client,
    document_map: DashMap<String, String>,
    tree_map: DashMap<String, Tree>,
    // TODO Symbols Map
}

use SemanticTokenType as ST;
pub const TOKEN_TYPES: [SemanticTokenType; 3] = [ST::KEYWORD, ST::TYPE, ST::FUNCTION];

#[tower_lsp::async_trait]
impl LanguageServer for Backend {
    async fn initialize(&self, _: InitializeParams) -> Result<InitializeResult> {
        // Ok(InitializeResult::default())
        Ok(InitializeResult {
            server_info: None,
            capabilities: ServerCapabilities {
                // inlay_hint_provider: Some(OneOf::Left(true)),
                text_document_sync: Some(TextDocumentSyncCapability::Kind(
                    TextDocumentSyncKind::FULL,
                )),
                // completion_provider: Some(CompletionOptions {
                //     resolve_provider: Some(false),
                //     trigger_characters: Some(vec![".".to_string()]), // TODO make ${ a trigger
                //     // character
                //     work_done_progress_options: Default::default(),
                //     all_commit_characters: None,
                //     completion_item: None,
                // }),
                // execute_command_provider: Some(ExecuteCommandOptions {
                //     commands: vec!["dummy.do_something".to_string()],
                //     work_done_progress_options: Default::default(),
                // }),
                workspace: Some(WorkspaceServerCapabilities {
                    workspace_folders: Some(WorkspaceFoldersServerCapabilities {
                        supported: Some(true),
                        change_notifications: Some(OneOf::Left(true)),
                    }),
                    file_operations: None,
                }),
                // TODO Consider adding more options like in the tower-lsp boilerplate example
                semantic_tokens_provider: Some(
                    SemanticTokensServerCapabilities::SemanticTokensOptions(
                        SemanticTokensOptions {
                            legend: SemanticTokensLegend {
                                token_types: TOKEN_TYPES.to_vec(),
                                token_modifiers: vec![],
                            },
                            ..Default::default()
                        },
                    ),
                ),
                // definition: Some(GotoCapability::default()),
                // definition_provider: Some(OneOf::Left(true)),
                // references_provider: Some(OneOf::Left(true)),
                // rename_provider: Some(OneOf::Left(true)),
                ..ServerCapabilities::default()
            },
        })
    }

    async fn initialized(&self, _: InitializedParams) {
        self.client
            .log_message(MessageType::INFO, "server initialized!")
            .await;
    }

    async fn shutdown(&self) -> Result<()> {
        Ok(())
    }

    async fn did_open(&self, params: DidOpenTextDocumentParams) {
        self.client
            .log_message(
                MessageType::INFO,
                format!("file opened: {}", params.text_document.uri),
            )
            .await;

        self.on_change(
            params.text_document.uri,
            params.text_document.text,
            params.text_document.version,
        )
        .await
    }
    async fn did_change(&self, mut params: DidChangeTextDocumentParams) {
        self.client
            .log_message(
                MessageType::INFO,
                format!("file changed: {}", params.text_document.uri),
            )
            .await;
        //
        // Not sure about this interface or method name, but it's in the template I'm hacking apart
        self.on_change(
            params.text_document.uri,
            std::mem::take(&mut params.content_changes[0].text),
            params.text_document.version,
        )
        .await
    }
}

// struct TextDocumentItem {
//     uri: Url,
//     text: String,
//     version: i32,
// }
impl Backend {
    // async fn on_change(&self, params: TextDocumentItem) {
    async fn on_change(&self, uri: Url, text: String, version: i32) {
        // example uses a rope, I'm just going to use a u8 vec or str.
        // let rope = ropey::Rope::from_str(&params.text);
        self.document_map.insert(uri.to_string(), text.clone());
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");

        let tree = parser.parse(&text, None).unwrap();
        // self.client
        //     .log_message(MessageType::INFO, format!("{:?}", errors))
        //     .await;

        // TODO combine this whole block into a single function
        let mut error_finder = ErrorFinder::new().unwrap();
        error_finder.find_syntax_errors(&tree, &text).unwrap();
        error_finder.find_missing_nodes(&tree).unwrap();
        let mut ident_finder = IdentiFinder::new(&tree, text.as_bytes()).unwrap();
        ident_finder.find_refs(&tree, text.as_bytes()).unwrap();
        ident_finder.find_defs(&tree, text.as_bytes()).unwrap();

        let mut diagnostics: Vec<Diagnostic> = error_finder
            .syntax_errors()
            .iter()
            .map(|e| e.clone().into())
            .collect();

        if let Err(v) = ident_finder.check_symbols() {
            diagnostics.extend(v.iter().map(|e| e.clone().into()))
        }
        self.client
            .publish_diagnostics(uri.clone(), diagnostics, Some(version))
            .await;

        // To stop a clippy warning about using format, and an async error
        // DEBUGGING
        // let sexp = &tree.root_node().to_sexp().to_string();
        // self.client.log_message(MessageType::INFO, sexp).await;

        // TODO pop old tree and update with new one
        self.tree_map.insert(uri.to_string(), tree);

        // self.client
        //     .log_message(MessageType::INFO, &format!("{:?}", semantic_tokens))
        //     .await;
        // self.semantic_token_map
        //     .insert(uri.to_string(), semantic_tokens);
    }
}

#[tokio::main]
async fn main() {
    let stdin = tokio::io::stdin();
    let stdout = tokio::io::stdout();

    let (service, socket) = LspService::new(|client| Backend {
        client,
        document_map: DashMap::new(),
        tree_map: DashMap::new(),
    });
    Server::new(stdin, stdout, socket).serve(service).await;
}
