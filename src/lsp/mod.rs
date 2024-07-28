//! Alternative implementation for the lsp using the `tower-lsp` crate.
use crate::ast::{
    ts_to_ast, Ast, CommandType, ComputeDef, FixDef, GenericCommand, NamedCommand, PartialAst,
};
use crate::check_commands;
use crate::check_styles::check_styles;
use crate::docs::docs_map::DOCS_MAP;
use crate::docs::DOCS_CONTENTS;
use crate::error_finder::ErrorFinder;
use crate::identifinder::{unused_variables, IdentiFinder};
use crate::lammps_errors::Warnings;
use crate::utils::get_symbol_at_point;
use dashmap::DashMap;
use std::str::FromStr;
use std::sync::Arc;
use tower_lsp::jsonrpc::Result;
use tower_lsp::lsp_types::*;
use tower_lsp::{Client, LanguageServer};
use tree_sitter::{Parser, Tree};

/// Core LSP Server Application
#[derive(Debug)]
pub struct Backend {
    client: Client,
    /// Map of file URI to document text
    document_map: DashMap<String, String>,
    /// Map of file URI to document tree
    tree_map: DashMap<String, Tree>,
    // TODO: Symbols Map
    /// Finds Symbols maps.
    /// Wrapped in RwLock for Interior Mutability
    identifinder: std::sync::RwLock<IdentiFinder>,

    ast: Arc<std::sync::RwLock<Ast>>,
}

impl Backend {
    pub fn new(client: Client) -> Self {
        Self {
            client,
            document_map: DashMap::new(),
            tree_map: DashMap::new(),
            identifinder: IdentiFinder::new_no_parse()
                .expect("Failed to compile Tree-sitter Query")
                .into(),
            ast: Arc::new(Ast::default().into()),
        }
    }
}

use SemanticTokenType as ST;
// TODO: Update tokens to have all the needed ones
// Match these to the highlight.scm in the tree-sitter package
pub const TOKEN_TYPES: [SemanticTokenType; 3] = [ST::KEYWORD, ST::TYPE, ST::FUNCTION];

#[tower_lsp::async_trait]
impl LanguageServer for Backend {
    async fn initialize(&self, _: InitializeParams) -> Result<InitializeResult> {
        Ok(InitializeResult {
            server_info: Some(ServerInfo {
                name: "LAMMPS Analyser".into(),
                version: None,
            }),
            capabilities: ServerCapabilities {
                text_document_sync: Some(TextDocumentSyncCapability::Kind(
                    // NOTE: Currently only supports full file update
                    TextDocumentSyncKind::FULL,
                )),
                workspace: Some(WorkspaceServerCapabilities {
                    workspace_folders: Some(WorkspaceFoldersServerCapabilities {
                        supported: Some(true),
                        change_notifications: Some(OneOf::Left(true)),
                    }),
                    file_operations: None,
                }),
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
                document_symbol_provider: Some(OneOf::Left(true)),
                workspace_symbol_provider: Some(OneOf::Left(true)),
                hover_provider: Some(HoverProviderCapability::Simple(true)),

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
        self.client
            .log_message(
                MessageType::INFO,
                format!("Tried to find symbols for {}", params.text_document.uri),
            )
            .await;

        // Assuming the identifier is properly updated
        // TODO: need to flatten???
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
            })
            .collect();

        Ok(Some(DocumentSymbolResponse::Flat(symbols)))
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
            .find(|(x, _)| *x == name_and_type);

        let Some(symbol) = symbol else {
            // For the case no identifier is found.
            // This should not happen, however, but this removes an unwrap.
            return Ok(None);
        };

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
        let point = params.text_document_position_params.position.into();

        // Unwrap should be ok. Panics if can't get lock result
        let ast = self.ast.read().unwrap();

        let command = ast.find_node(point);

        // FIXME: For unknown fix/compute styles does some weird stuff.
        if let Some(command) = command {
            let doc_name = match &command.command_type {
                CommandType::NamedCommand(NamedCommand::Fix(FixDef { fix_style, .. })) => {
                    DOCS_MAP.fixes().get(fix_style)
                }
                CommandType::NamedCommand(NamedCommand::Compute(ComputeDef {
                    compute_style,
                    ..
                })) => DOCS_MAP.computes().get(compute_style),

                CommandType::GenericCommand(GenericCommand { name, .. }) => {
                    DOCS_MAP.commands().get(name)
                }

                _ => None,
            };

            let Some(doc_name) = doc_name else {
                return Ok(None);
            };

            // Read the actual docs
            let hover_text = DOCS_CONTENTS
                .get(format!("{doc_name}.md").as_str())
                .expect("No docs for this command")
                .to_string();

            Ok(Some(Hover {
                contents: HoverContents::Markup(MarkupContent {
                    kind: MarkupKind::Markdown,
                    value: hover_text,
                }),
                range: None,
            }))
        } else {
            Ok(None)
        }
    }
}

impl Backend {
    // async fn on_change(&self, params: TextDocumentItem) {
    /// Main method that is called on change to the file.
    ///
    /// Currently the only method used for updating the server state.
    ///
    /// Actions performed:
    /// - Parse file with tree-sitter
    /// - Convert the TS tree into a `Ast`
    /// - Publish errors in the `Ast`a
    /// - Find any identifiers in the file.
    /// - Run checks on the fixes and publish these errors.
    /// - Checks for any issues with the symbols of the file.
    async fn on_change(&self, uri: Url, text: String, version: i32) {
        self.document_map.insert(uri.to_string(), text);

        let text = self.document_map.get(&uri.to_string()).unwrap();

        // TODO: Store parser in the state of the backend.
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");

        let tree = parser.parse(&*text, None).unwrap();

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

        *self.ast.write().unwrap() = ast;

        // TODO: combine this whole block into a single function
        let mut error_finder = ErrorFinder::new().unwrap();
        error_finder.find_syntax_errors(&tree, &*text).unwrap();
        error_finder.find_missing_nodes(&tree).unwrap();
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

        let mut diagnostics: Vec<Diagnostic> = error_finder
            // TODO: implement `From<&SyntaxError>` for `Diagnostic` so explicit cloning isn't needed.
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

        // TODO: pop old tree and update with new one
        self.tree_map.insert(uri.to_string(), tree);
    }
}
