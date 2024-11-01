use std::fmt;
use crate::topology::*;

impl fmt::Display for Simplex{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
     	write!(f, "{}-simplex '{:?}'", self.dimension, self.value)
    }
}
impl fmt::Display for SimplicialComplex{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    	write!(f, "simplicial {}-complex {{\n", self.dimension);
    	for simplex in &self.simplicies{
    		write!(f, "\t{}\n", simplex);
    	}
     	write!(f, "}}")
    }
}

impl fmt::Display for Ring{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
     	match self{
     		Ring::R => write!(f, "R"),
     		Ring::Z => write!(f, "Z"),
     		Ring::TRIVIAL => write!(f,"0")
     	}
    }
}
