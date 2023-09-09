use std::io::stdin;

use rust_systre::spacegroups::parse_table::parse_space_group_table;


fn main() {
    let tables = parse_space_group_table(stdin()).unwrap();
    let settings = tables.settings;

    println!("{}", settings.get("R3:R").unwrap());
    println!("{}", settings.get("P4132").unwrap());
}
