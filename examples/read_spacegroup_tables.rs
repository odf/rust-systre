use std::io::stdin;

use rust_systre::spacegroups::parse_table::parse_space_group_table;


fn main() {
    let tables = parse_space_group_table(stdin()).unwrap();
    let settings = tables.settings;
    let alias = tables.alias;
    let lookup = tables.lookup;

    print!("{}\n\n", settings.get("c2mm").unwrap());
    print!("{}\n\n", settings.get("R3:R").unwrap());
    print!("{}\n\n", settings.get("P4132").unwrap());

    let not_found: String = "<not found>".into();
    for name in ["Pm3", "A2/a", "XYZ"] {
        println!("{} -> {}", name, alias.get(name).unwrap_or(&not_found));
    }
    println!();

    println!();
    for i in [10, 20, 30] {
        print!("{}\n\n", lookup[i]);
    }
}
