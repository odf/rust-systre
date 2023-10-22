use rust_systre::spacegroups::parse_table::{parse_space_group_table, Tables};

#[macro_use]
extern crate lazy_static;

lazy_static! {
    static ref TABLES: Tables = parse_space_group_table(
        include_str!(
            concat!(env!("CARGO_MANIFEST_DIR"), "/src/data/sgtable.data")
        ).as_bytes()
    ).unwrap();
}


fn main() {
    let settings = &TABLES.settings;
    let alias = &TABLES.alias;
    let lookup = &TABLES.lookup;

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
