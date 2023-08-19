#[derive(Debug, PartialEq)]
pub enum CrystalSystem2d {
    Oblique,
    Rectangular,
    Square,
    Hexagonal,
}

pub enum CrystalSystem3d {
    Cubic,
    Orthorhombic,
    Hexagonal,
    Tetragonal,
    Trigonal,
    Monoclinic,
    Triclinic,
}

pub enum Centering2d {
    Primitive,
    Centered,
}

pub enum Centering3d {
    Primitive,
    FaceCentered,
    BodyCentered,
    Rhombohedral,
    AFaceCentered,
    BFaceCentered,
    CFaceCentered,
}

pub struct SpaceGroup2d {
    dimension: usize,
    crystal_system: CrystalSystem2d,
    centering: Centering2d,
    full_name: String,
    group_name: String,
    extension: String,
}
