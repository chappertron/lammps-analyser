macro_rules! derive_styles {
    // Comma inside repeat requires trailing comma
    ($EnumName:ident => $(($variant:tt, $lit:literal),)+) => {

        // #[allow(non_camel_case_types)]
        #[derive(Debug,Eq,PartialEq,Clone,Hash)]
        #[non_exhaustive]
        pub enum $EnumName {
            $($variant,)+
            InvalidStyle,
        }


        impl From<&str> for $EnumName {
            fn from(value: &str) -> Self {
                match value {
                    $($lit => $EnumName::$variant,)+
                    _ => $EnumName::InvalidStyle,
                }
            }
        }
        impl std::fmt::Display for $EnumName {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match self {
                    $($EnumName::$variant => write!(f, $lit),)+
                    $EnumName::InvalidStyle => write!(f, "InvalidStyle"),
                }
            }
        }

        impl Default for $EnumName {
            fn default() -> Self {
                $EnumName::InvalidStyle
            }
        }
    };
}

pub mod pair_styles;

// derive_styles!(
//     PairStyles =>
//     (Style1, "style1"),
//     (Style2, "style2"),
//     (Style3, "style3"),
//     (Style4, "style4"),
// );
