flamegraph:
	CARGO_PROFILE_RELEASE_DEBUG=true cargo flamegraph -b lammps-analyser --root -- ./example_input_scripts/in.nemd
typos:
	typos --exclude 'lammps_docs*/' --exclude 'example_input_scripts' --format brief

