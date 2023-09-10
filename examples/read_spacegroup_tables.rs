use std::io::stdin;

use rust_systre::spacegroups::parse_table::parse_space_group_table;


fn main() {
    let tables = parse_space_group_table(stdin()).unwrap();
    let settings = tables.settings;

    print!("{}\n\n", settings.get("c2mm").unwrap());
    print!("{}\n\n", settings.get("R3:R").unwrap());
    print!("{}\n\n", settings.get("P4132").unwrap());
}
