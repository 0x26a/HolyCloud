use crate::topology::*;
use nalgebra_lapack::SVD;

pub struct Homology<const T: usize>{
	pub start: f64,
	pub end: f64,
	pub ranks: [Rank;T],
	pub torsions: [Torsion;T]
}
impl SimplicialComplex{
	// A_module H_n(X,A)
	pub fn H(&self, n: usize) -> (usize, Torsion){
		if n <= self.dimension{
			return match self.ring{
				Ring::R => {
					let hom_rank = self.b(n);
					if hom_rank > 0 { (hom_rank, vec![]) } else { (0, vec![]) }
				},
				Ring::Z => {
					let mut torsion: Vec<usize> = Vec::new();
					let hom_rank = self.b(n);
					for invariant in smith_normal_form::compute(self.boundary(n+1)){
						if invariant >= 2.0 { torsion.push(invariant as usize) }
					}
					(hom_rank, torsion)

				},
				_ => panic!("H(-,{}) homology isn't supported", self.ring)
			};
		}
		(0, vec![])
	}

	// Betti number b_n(X)
	pub fn b(&self,n:usize) -> usize{
		if n <= self.dimension {
			self.dim_C(n) - SVD::new(self.boundary(n)).unwrap().rank(TOL) - SVD::new(self.boundary(n+1)).unwrap().rank(TOL) 
		}
		else { 0 }
	}
}


