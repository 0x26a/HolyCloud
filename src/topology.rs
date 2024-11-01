pub mod samples;
mod display;
pub mod homology;
mod smith_normal_form;

use nalgebra::DMatrix;

const TOL: f64 = 0.000001;

#[derive(Clone)]
#[derive(Debug)]
#[derive(PartialEq)]
pub enum Ring{
	TRIVIAL,
	R, 
	Z, 
}
pub type Rank = usize;
pub type Torsion = Vec<usize>;

#[derive(Clone)]
#[derive(Debug)]
pub struct const_Simplex<'a>{
	pub value: &'a [usize],
	pub dimension: usize
}
#[derive(Clone)]
#[derive(Debug)]
pub struct Simplex{
	pub value: Vec<usize>,
	pub dimension: usize
}

#[derive(Debug)]
pub struct SimplicialComplex{
	pub ring: Ring,
	simplicies : Vec<Simplex>,
	pub dimension: usize,
	chains: Vec<Vec<Simplex>>

}

impl Simplex{
	pub fn from(value: Vec<usize>) -> Simplex{ let l = value.len(); Simplex{ value: value, dimension: l - 1 } }
}

impl SimplicialComplex{
	pub fn from(simplicies: Vec<Simplex>, ring: Ring) -> SimplicialComplex{ 
		let mut s = SimplicialComplex{ ring, simplicies, dimension: 0, chains: vec![] };
		s.gen_C_spaces();
		return s;
	}	
	pub fn const_from(const_simplicies: &[&const_Simplex], ring: Ring) -> SimplicialComplex{
		let mut simplicies: Vec<Simplex> = Vec::new();
		for simp in const_simplicies{ simplicies.push(Simplex{ value: simp.value.to_vec(), dimension: simp.dimension}); }
		let mut s = SimplicialComplex{ ring, simplicies, dimension: 0, chains: vec![] };
		s.gen_C_spaces();
		return s;
	}

	fn gen_C_spaces(&mut self){
		let mut l = 0;
		let mut dimension = 0;
		let mut chains: Vec<Vec<Simplex>> = vec![vec![]];
		for f in &self.simplicies{
			while l < f.dimension{
				chains.push(vec![]);
				l += 1;
			}
			chains[f.dimension].push(f.clone());
			dimension = f.dimension;
		}
		self.dimension = dimension;
		self.chains = chains;
	}


	pub fn dim_C(&self, n: usize) -> usize{ self.chains[n].len() }

	pub fn boundary(&self, n: usize) -> DMatrix<f64>{
		match n{
			
			n if n == 0 => DMatrix::from_vec(1,self.chains[0].len(),vec![0.0;self.chains[0].len()]),

			n if n > self.dimension => DMatrix::from_vec(1,1,vec![0.0]),

			_ => {
				let p = self.dim_C(n);
				let q = self.dim_C(n-1);
				let mut vec_of_mat: Vec<f64> = Vec::new();
				for simplex in self.chains[n].iter(){
					let mut coord = self.boundary_basis(n, q, simplex);
					vec_of_mat.append(&mut coord);
				}
				DMatrix::from_vec(q,p,vec_of_mat)
			}
		}
	}
	fn boundary_basis(&self, n: usize, p: usize, simplex: &Simplex) -> Vec<f64>{

		let mut position = 0;
		let mut vector = vec![0.0;p];
		let mut value = simplex.value.clone();
		let mut i = 0;
		let l = value.len();

		while i < l{
			value.remove(i);
			for s in &self.chains[n-1]{
				if s.value == value{ break; }
				position += 1
			}
			vector[position] = if i % 2 == 0 { 1.0 } else { -1.0 };
			value = simplex.value.clone();
			position = 0;
			i += 1;
		}
		vector
	}


}

