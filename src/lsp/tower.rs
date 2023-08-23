use dashmap::DashMap;
use lammps_analyser::error_finder::{self, ErrorFinder};
///! Alternative implementation for the lsp using the `tower-lsp` crate.
use tower_lsp::jsonrpc::Result;
use tower_lsp::lsp_types::*;
use tower_lsp::{Client, LanguageServer, LspService, Server};
use tree_sitter::{Parser, Tree};
use tree_sitter_lammps;
#[derive(Debug)]
struct Backend {
    client: Client,
    document_map: DashMap<String, String>,
    tree_map: DashMap<String, Tree>,
}

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
                semantic_tokens_provider: None,
                // definition: Some(GotoCapability::default()),
                definition_provider: Some(OneOf::Left(true)),
                references_provider: Some(OneOf::Left(true)),
                rename_provider: Some(OneOf::Left(true)),
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

        // Not sure about this interface or method name, but it's in the template I'm hacking apart
        self.on_change(TextDocumentItem {
            uri: params.text_document.uri,
            text: params.text_document.text,
            version: params.text_document.version,
        })
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
        self.on_change(TextDocumentItem {
            uri: params.text_document.uri,
            text: std::mem::take(&mut params.content_changes[0].text),
            version: params.text_document.version,
        })
        .await
    }
}
struct TextDocumentItem {
    uri: Url,
    text: String,
    version: i32,
}
impl Backend {
    async fn on_change(&self, params: TextDocumentItem) {
        // example uses a rope, I'm just going to use a u8 vec or str.
        // let rope = ropey::Rope::from_str(&params.text);
        let text = &params.text;
        self.document_map
            .insert(params.uri.to_string(), text.clone());
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");

        let tree = parser.parse(text, None).unwrap();
        // self.client
        //     .log_message(MessageType::INFO, format!("{:?}", errors))
        //     .await;
        let mut error_finder = ErrorFinder::new().unwrap();
        error_finder.find_syntax_errors(&tree, &text).unwrap();
        error_finder.find_missing_nodes(&tree).unwrap();
        // COMPILES BUT INCORRECT BEHAVIOUR
        let diagnostics = error_finder
            .syntax_errors()
            .iter()
            .map(|e| e.clone().into())
            .collect();

        self.client
            .publish_diagnostics(params.uri.clone(), diagnostics, Some(params.version))
            .await;

        self.client
            .log_message(
                MessageType::INFO,
                &format!("{}", tree.root_node().to_sexp()),
            )
            .await;
        // TODO pop old tree and update with new one
        self.tree_map.insert(params.uri.to_string(), tree);

        // self.client
        //     .log_message(MessageType::INFO, &format!("{:?}", semantic_tokens))
        //     .await;
        // self.semantic_token_map
        //     .insert(params.uri.to_string(), semantic_tokens);
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
