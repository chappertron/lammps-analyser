use std::collections::HashMap;

use once_cell::sync::Lazy;

pub mod docs_map;

// Generated in the build script
const ALL_DOCS: &[(&str, &str)] = &include!(concat!(env!("OUT_DIR"), "/all_the_docs.rs"));

/// Map between the docs file name and the actual contents.
/// Baked into the binary to avoid folder lookup issues.
pub static DOCS_CONTENTS: Lazy<HashMap<&str, &str>> =
    Lazy::new(|| HashMap::from_iter(ALL_DOCS.iter().map(|&(a, b)| (a, b))));
