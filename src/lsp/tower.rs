//! Alternative implementation for the lsp using the `tower-lsp` crate.
use dashmap::DashMap;
use lammps_analyser::ast::{ts_to_ast, Ast, CommandType, NamedCommand, PartialAst};
use lammps_analyser::check_commands;
use lammps_analyser::check_styles::check_styles;
use lammps_analyser::error_finder::ErrorFinder;
use lammps_analyser::identifinder::{unused_variables, IdentiFinder};
use lammps_analyser::lammps_errors::Warnings;
use lammps_analyser::utils::{get_symbol_at_point, position_to_point};
use std::str::FromStr;
use std::sync::Arc;
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
    /// Map of uri to document text
    document_map: DashMap<String, String>,
    /// Map of uri to document tree
    tree_map: DashMap<String, Tree>,
    // TODO: Symbols Map
    /// Finds Symbols maps.
    /// Wrapped in RwLock for Interior Mutability
    identifinder: std::sync::RwLock<IdentiFinder>,

    ast: Arc<std::sync::RwLock<Option<Ast>>>,
}

use SemanticTokenType as ST;
// TODO: Update tokens to have all the needed ones
// Match these to the highlight.scm in the tree-sitter package
pub const TOKEN_TYPES: [SemanticTokenType; 3] = [ST::KEYWORD, ST::TYPE, ST::FUNCTION];

#[tower_lsp::async_trait]
impl LanguageServer for Backend {
    async fn initialize(&self, _: InitializeParams) -> Result<InitializeResult> {
        // Ok(InitializeResult::default())
        Ok(InitializeResult {
            server_info: Some(ServerInfo {
                name: "LAMMPS Analyser".into(),
                version: None,
            }),
            capabilities: ServerCapabilities {
                // inlay_hint_provider: Some(OneOf::Left(true)),
                text_document_sync: Some(TextDocumentSyncCapability::Kind(
                    // Currently only supports full file update
                    TextDocumentSyncKind::FULL,
                )),
                // completion_provider: Some(CompletionOptions {
                //     resolve_provider: Some(false),
                //     trigger_characters: Some(vec![".".to_string()]), // TODO: make ${ a trigger
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
                // TODO: Consider adding more options like in the tower-lsp boilerplate example
                semantic_tokens_provider: Some(
                    SemanticTokensServerCapabilities::SemanticTokensOptions(
                        SemanticTokensOptions {
                            legend: SemanticTokensLegend {
                                token_types: TOKEN_TYPES.into(),
                                token_modifiers: vec![],
                            },
                            ..Default::default()
                        },
                    ),
                ),
                definition_provider: Some(OneOf::Left(true)),
                // TODO: Add options?
                document_symbol_provider: Some(OneOf::Left(true)),
                workspace_symbol_provider: Some(OneOf::Left(true)),
                hover_provider: Some(HoverProviderCapability::Simple(true)),

                // definition: Some(GotoCapability::default()),
                // definition_provider: Some(OneOf::Left(true)),
                // references_provider: Some(OneOf::Left(true)),
                // rename_provider: Some(OneOf::Left(true)),
                ..Default::default()
            },
        })
    }

    async fn initialized(&self, _: InitializedParams) {
        self.client
            .log_message(MessageType::INFO, "Server Initialized!")
            .await;
    }

    async fn shutdown(&self) -> Result<()> {
        Ok(())
    }

    async fn did_open(&self, params: DidOpenTextDocumentParams) {
        self.client
            .log_message(
                MessageType::INFO,
                format!("File opened: {}", params.text_document.uri),
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
                format!("File changed: {}", params.text_document.uri),
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

    async fn symbol(
        &self,
        params: WorkspaceSymbolParams,
    ) -> Result<Option<Vec<SymbolInformation>>> {
        // TODO: This should use the query for the document_symbols to find the workspace symbols
        // TODO: cache the symbols rather than refinding them all?
        let mut symbols = vec![];
        for r in &self.document_map {
            let uri = r.key();
            let params = DocumentSymbolParams {
                text_document: TextDocumentIdentifier {
                    uri: Url::from_str(uri).expect("Invalid uri"),
                },
                work_done_progress_params: params.work_done_progress_params.clone(),
                partial_result_params: params.partial_result_params.clone(),
            };

            symbols.extend(
                self.document_symbol(params)
                    .await?
                    .into_iter()
                    .flat_map(|x| match x {
                        DocumentSymbolResponse::Flat(x) => x,
                        DocumentSymbolResponse::Nested(_) => {
                            panic!("Nested symbols not supported???")
                        }
                    }),
            );
        }

        Ok(Some(symbols))
    }
    async fn document_symbol(
        &self,
        params: DocumentSymbolParams,
    ) -> Result<Option<DocumentSymbolResponse>> {
        // Ok(Some(DocumentSymbolResponse::Flat(vec![])))
        self.client
            .log_message(
                MessageType::ERROR,
                format!("Tried to find symbols for {}", params.text_document.uri),
            )
            .await;

        // Assuming the identifier is properly updated
        // let symbols = ;
        // need to flatten???
        #[allow(deprecated)] // Used for the `deprecated` field of the `SymbolInformation` struct
        let symbols = self
            .identifinder
            .read()
            .unwrap()
            .symbols()
            .iter()
            // TODO: properly match this field to the compute/var/fix type
            // TODO: Set correct Positions for the Symbols
            .flat_map(|(x, v)| {
                v.defs().iter().map(|s| SymbolInformation {
                    name: x.name.clone(),
                    kind: SymbolKind::from(s.ident_type),
                    location: Location {
                        uri: params.text_document.uri.clone(),
                        range: s.range().into_lsp_types(),
                    },
                    container_name: None,
                    tags: None,

                    deprecated: None,
                })
                // .collect::<Vec<_>>()
            })
            .collect();

        // self.client
        //     .log_message(
        //         MessageType::INFO,
        //         format!(
        //             "Found Symbols for {}: {:?}",
        //             params.text_document.uri, &symbols
        //         ),
        //     )
        //     .await;

        Ok(Some(DocumentSymbolResponse::Flat(symbols)))
        // todo!()
    }

    async fn goto_definition(
        &self,
        params: GotoDefinitionParams,
    ) -> Result<Option<GotoDefinitionResponse>> {
        let uri = params.text_document_position_params.text_document.uri;
        let position = params.text_document_position_params.position;
        let point = position.into();

        self.client
            .log_message(MessageType::INFO, format!("At point: {point:?}"))
            .await;

        let identifiers = self.identifinder.read().unwrap();

        let Some(name_and_type) = get_symbol_at_point(&identifiers, &point) else {
            return Ok(None);
        };

        let symbol = identifiers
            .symbols()
            .iter()
            .find(|(x, _)| *x == name_and_type)
            .expect("Could not find symbol"); // TODO: Bubble this up

        Ok(Some(GotoDefinitionResponse::Array(
            symbol
                .1
                .defs()
                .iter()
                .map(|def| Location {
                    uri: uri.clone(),
                    range: def.range().into_lsp_types(),
                })
                .collect(),
        )))
    }

    async fn hover(&self, params: HoverParams) -> Result<Option<Hover>> {
        let uri = params.text_document_position_params.text_document.uri;
        let ts_point = position_to_point(&params.text_document_position_params.position);

        let tree = self.tree_map.get(&uri.to_string()).unwrap();
        let mut tree_cursor = tree.walk();

        let point = params.text_document_position_params.position.into();

        // TODO: Tidy the unwraps
        let ast = self.ast.read().unwrap();

        let command = ast
            .as_ref()
            .map(|ast| ast.find_node(point))
            .flatten()
            .unwrap();

        let hover_text = format!(
            "**Hovering** over: {:?}\n**Raw Position: {:?}**",
            command, point
        );
        // Old hard coded example. TODO: remove
        // include_str!("../../lammps_docs_md/fix_langevin.md").to_owned()
        Ok(Some(Hover {
            contents: HoverContents::Markup(MarkupContent {
                kind: MarkupKind::Markdown,
                // TODO: Actually find the right document for the command.
                value: hover_text,
            }),
            range: None,
        }))
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
        self.document_map.insert(uri.to_string(), text);

        let text = self.document_map.get(&uri.to_string()).unwrap();

        // TODO: Store parser in the state
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");

        let tree = parser.parse(&*text, None).unwrap();
        // self.client
        //     .log_message(MessageType::INFO, format!("{:?}", errors))
        //     .await;

        let ast = ts_to_ast(&tree, text.as_bytes());

        let ast = match ast {
            Ok(ast) => ast,
            Err(PartialAst { ast, errors }) => {
                for e in errors {
                    println!("{}", e);

                    self.client
                        .publish_diagnostics(
                            uri.clone(),
                            vec![Diagnostic {
                                severity: Some(DiagnosticSeverity::ERROR),
                                message: format!("Error Parsing TS To AST: {}", e),
                                range: e.span.into_lsp_types(),
                                ..Default::default()
                            }],
                            None,
                        )
                        .await;
                }

                ast
            }
        };

        *self.ast.write().unwrap() = Some(ast);

        // TODO: combine this whole block into a single function
        let mut error_finder = ErrorFinder::new().unwrap();
        error_finder.find_syntax_errors(&tree, &*text).unwrap();
        error_finder.find_missing_nodes(&tree).unwrap();
        // let mut ident_finder = IdentiFinder::new(&tree, text.as_bytes()).unwrap();
        // ident_finder.find_refs(&tree, text.as_bytes()).unwrap();
        // ident_finder.find_defs(&tree, text.as_bytes()).unwrap();
        self.identifinder
            .write()
            .unwrap()
            .find_symbols(&tree, text.as_bytes())
            .unwrap();

        let invalid_styles = check_styles(&tree, text.as_bytes());

        // Parsing fixes

        let fix_errors = self
            .ast
            .read()
            .expect("Taking read lock on ast")
            .as_ref()
            .unwrap()
            .commands
            .iter()
            .filter_map(|command| {
                if let CommandType::NamedCommand(NamedCommand::Fix(fix)) = &command.command_type {
                    Some(check_commands::fixes::check_fix(fix))

                // TODO: add checking for computes here.
                } else {
                    None
                }
            })
            .filter_map(|x| x.err())
            .collect::<Vec<_>>();

        dbg!(&fix_errors);

        let mut diagnostics: Vec<Diagnostic> = error_finder
            // TODO: add methods to create ownership, so clowning isn't needed???
            // TODO: implement into for vec of error
            .syntax_errors()
            .iter()
            .map(|e| e.clone().into())
            .collect();

        if let Err(v) = self.identifinder.read().unwrap().check_symbols() {
            diagnostics.extend(v.into_iter().map(|e| e.into()))
        }
        if let Ok(v) = invalid_styles {
            diagnostics.extend(v.into_iter().map(|e| e.into()))
        }

        // Add unused symbols
        diagnostics.extend(
            unused_variables(self.identifinder.read().unwrap().symbols())
                .into_iter()
                .map(|x| Warnings::from(x).into()),
        );

        // Zero or more
        diagnostics.extend(fix_errors.into_iter().map(|e| e.into()));

        self.client
            .publish_diagnostics(uri.clone(), diagnostics, Some(version))
            .await;

        // TODO: Could be excessive
        // TODO: handle result
        // self.client.workspace_diagnostic_refresh().await;

        // DEBUGGING
        // To stop a clippy warning about using format, and an async error
        // let sexp = &tree.root_node().to_sexp().to_string();
        // self.client.log_message(MessageType::INFO, sexp).await;

        // TODO: pop old tree and update with new one
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
        identifinder: IdentiFinder::new_no_parse().unwrap().into(),
        ast: Arc::new(None.into()),
    });
    Server::new(stdin, stdout, socket).serve(service).await;
}
