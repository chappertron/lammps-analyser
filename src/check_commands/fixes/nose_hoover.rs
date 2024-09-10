use crate::{
    ast::{ArgumentKind, FixDef},
    check_commands::{invalid_arguments, utils},
    fix_styles::FixStyle,
};

/// Generic Parsing of the Nose-Hoover Fixes
/// TODO: Finish ME
/// TODO: Error if a temp or barostat keyword is not being used
pub(crate) fn parse_nh_fixes(fix: &FixDef) -> Result<(), invalid_arguments::InvalidArgumentsType> {
    let args = &fix.args;

    // Iterate through and check validity of the argument
    let mut iter = args.iter();

    // TODO: Add the more niche fixes to these.
    let is_barostatting = matches!(fix.fix_style, FixStyle::Npt | FixStyle::Nph);
    let is_thermostatting = matches!(fix.fix_style, FixStyle::Npt | FixStyle::Nvt);
    // Try either the LAMMPS way or to use peek with a while let loop?

    // Errors if the keyword is not being used with  barostatting style
    let barostat_only = |kwarg: &str| {
        if !is_barostatting {
            Err(invalid_arguments::InvalidArgumentsType::InvalidKeyword {
                kwarg: kwarg.to_string(),
                fix_style: fix.fix_style,
            })
        } else {
            Ok(())
        }
    };

    // Errors if the keyword is not being used with  barostatting style
    let thermostat_only = |kwarg: &str| {
        if !is_thermostatting {
            Err(invalid_arguments::InvalidArgumentsType::InvalidKeyword {
                kwarg: kwarg.to_string(),
                fix_style: fix.fix_style,
            })
        } else {
            Ok(())
        }
    };

    while let Some(arg) = iter.next() {
        match &arg.kind {
            // TODO: Change ArgName again if this changes again
            ArgumentKind::Word(kwarg) if kwarg == "temp" => {
                // TODO: check if there are 3 more elements
                thermostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 3, "<Tstart> <Tstop> <Tdamp>")?;
            }
            ArgumentKind::Word(kwarg)
                if matches!(
                    kwarg.as_ref(),
                    "iso" | "aniso" | "tri" | "x" | "y" | "z" | "xy" | "xz" | "yz"
                ) =>
            {
                barostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 3, "<Pstart> <Pstop> <Pdamp>")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "couple" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_enum(
                    &mut iter,
                    kwarg,
                    1,
                    &["none", "xyz", "xy", "xz", "yz"],
                )?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "tchain" => {
                thermostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 1, "<N> chain.")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "pchain" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 1, "<N> chain.")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "mtk" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "tloop" => {
                utils::kwarg_expected_floats(&mut iter, kwarg, 1, "<N> sub-cycles")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "ploop" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 1, "<N> sub-cycles")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "nreset" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 1, "<N> reset")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "drag" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 1, "<Df> drag factor")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "ptemp" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 1, "<Ttarget>")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "dilate" => {
                barostat_only(kwarg)?;
                // TODO: convert to a group
                utils::kwarg_expected_str(&mut iter, kwarg, 1, "<dilate-group-ID>")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "scalexy" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "scaleyz" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "scalexz" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }

            ArgumentKind::Word(kwarg) if kwarg == "flip" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }

            ArgumentKind::Word(kwarg) if kwarg == "fixedpoint" => {
                barostat_only(kwarg)?;
                utils::kwarg_expected_floats(&mut iter, kwarg, 1, "<x> <y> <z>")?;
            }
            ArgumentKind::Word(kwarg) if kwarg == "update" => {
                utils::kwarg_expected_enum(&mut iter, kwarg, 1, &["dipole", "dipole/dlm"])?
            }

            ArgumentKind::Word(kwarg) => Err(invalid_arguments::InvalidArgumentsType::Custom(
                format!("Unknown kwarg argument: {kwarg}",),
            ))?,

            _ => Err(invalid_arguments::InvalidArgumentsType::Custom(format!(
                "Unknown argument: {arg}",
            )))?,
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ast::from_node::FromNode, utils::testing::parse};

    #[test]
    fn valid_nvt() {
        let text = "fix NVT all nvt temp 1 $(v_T*1.5) $(100*dt)";
        let tree = parse(text);
        let node = tree.root_node().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert_eq!(fix.fix_id.name, "NVT");
        assert_eq!(fix.group_id.contents, "all");
        assert_eq!(fix.fix_style, FixStyle::Nvt);
        assert!(!fix.args.is_empty());

        assert_eq!(parse_nh_fixes(&fix), Ok(()));
    }

    #[test]
    fn valid_npt() {
        let text = "fix NPT all npt temp 1 $(v_T*1.5) $(100*dt) iso 1 $(v_p*1.5) ${pdamp}";
        let tree = parse(text);
        let node = tree.root_node().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert_eq!(fix.fix_id.name, "NPT");
        assert_eq!(fix.group_id.contents, "all");
        assert_eq!(fix.fix_style, FixStyle::Npt);
        assert!(!fix.args.is_empty());

        assert_eq!(dbg!(parse_nh_fixes(&fix)), Ok(()));
    }

    #[test]
    fn complicated_nph() {
        let text = "fix 2 ice nph x 1.0 1.0 0.5 y 2.0 2.0 0.5 z 3.0 3.0 0.5 yz 0.1 0.1 0.5 xz 0.2 0.2 0.5 xy 0.3 0.3 0.5 nreset 1000";
        let tree = parse(text);
        let node = tree.root_node().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert_eq!(fix.fix_id.name, "2");
        assert_eq!(fix.group_id.contents, "ice");
        assert_eq!(fix.fix_style, FixStyle::Nph);
        assert!(!fix.args.is_empty());

        assert_eq!(dbg!(parse_nh_fixes(&fix)), Ok(()));
    }

    #[test]
    fn invalid_nvt() {
        let text = "fix NVT all nvt temp 1.0 $(v_T*1.5) TEMP";
        let tree = parse(text);
        let node = tree.root_node().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert_eq!(fix.fix_id.name, "NVT");
        assert_eq!(fix.group_id.contents, "all");
        assert_eq!(fix.fix_style, FixStyle::Nvt);
        assert!(!fix.args.is_empty());

        // TODO: Turn into an error...
        assert!(parse_nh_fixes(&fix).is_err());
    }
}
