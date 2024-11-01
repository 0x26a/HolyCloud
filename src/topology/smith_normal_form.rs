use nalgebra::DMatrix;
use crate::topology::smith_normal_form::SNF_STATE::*;

pub enum SNF_STATE{
	NextStep(f64),
	GoToInitial,
	Halt(&'static str)
}
pub struct SmithNormalForm<'a>{
	m: &'a mut DMatrix<f64>,
	p: usize,
	q: usize
}
impl SmithNormalForm<'_>{
	pub fn block_diagonal(&mut self, a: f64){
		for i in 0..self.p{ if i > 0{ self.e1_ij(i,0,(-1.0)*self.m[(i,0)]/a); } }
		for j in 0..self.q{ if j > 0{ self.e2_ij(j,0,(-1.0)*self.m[(0,j)]/a); } }
	}
	pub fn bad_entry(&self, a: f64) -> Option<(usize,usize)>{
		for j in 0..self.q{
			if self.m[(0,j)].rem_euclid(a) != 0.0{
				return Some((0,j))
			}
		}

		for i in 0..self.p{
			if self.m[(i,0)].rem_euclid(a) != 0.0{
				return Some((i,0))
			}
		}
		for i in 0..self.p{
			for j in 0..self.q{
				if self.m[(i,j)].rem_euclid(a) != 0.0{
					return Some((i,j))
				}
			}
		}
		None
	}
	pub fn swap_to_top(&mut self, (i,j): (usize, usize)){
		self.m.swap_rows(i,0);
		self.m.swap_columns(j,0);
	}
	pub fn norm_min_pos(&self) -> Option<((usize, usize), f64)>{
		let mut i = 0;
		let mut j = 0;
		let mut x = self.m[(0,0)].abs();

		let mut x_no_abs = self.m[(0,0)];

		for _i in 0..self.p{
			for _j in 0..self.q{
				let v = self.m[(_i,_j)].abs();
				if (v < x && v > 0.0) || (x == 0.0){
					i = _i;
					j = _j;
					x = v;
					x_no_abs = self.m[(_i,_j)];
				}
			}
		}
		if x == 0.0{ return None; }
		Some(((i,j),x_no_abs))
	}
	pub fn e1_ij(&mut self, i: usize, j:usize, r: f64){
		self.m.set_row(i,&(self.m.row(i)+self.m.row(j)*r));
	}
	pub fn e2_ij(&mut self, i: usize, j:usize, r:f64){
		self.m.set_column(i,&(self.m.column(i) + self.m.column(j)*r));
	
	}

}

pub fn step(snf: &mut SmithNormalForm) -> SNF_STATE{
	//println!("init: {}", snf.m);

	match snf.norm_min_pos(){
		None => Halt("halted: M=0"),
		Some((pos,a)) => {
			snf.swap_to_top(pos);
			match snf.bad_entry(a){
				None => {
					snf.block_diagonal(a);
					NextStep(a)
				},
				Some((0,j)) => {
					let q = snf.m[(0,j)].div_euclid(a);
					snf.e2_ij(j,0,-q);
					GoToInitial
				},
				Some((i,0)) => {
					let q = snf.m[(i,0)].div_euclid(a);
					////println!("{}={}*{} + {}", snf.m[(i,0)],q,a,r);
					snf.e1_ij(i,0,-q);
					GoToInitial
				},
				Some((i,j)) => {
					snf.block_diagonal(a);
					match snf.bad_entry(a){
						None => NextStep(a),
						_ => {
							snf.e1_ij(0,i,1.0);
							NextStep(a)
						}
					}
				}
			}
		}
	}
}


pub fn compute(mut m: DMatrix<f64>) -> Vec<f64>{
	let mut invariants: Vec<f64> = Vec::new();
	let mut q = m.row(0).len();
	let mut p = m.column(0).len();
	if q == 1 && p == 1 && m[(0,0)] == 0.0{
		return vec![];
	}

	while p > 0 && q > 0 {
		match step(&mut SmithNormalForm{
			m: &mut m,
			p: p,
			q: q
		}){
			NextStep(a) => {
				invariants.push(a.abs());
				m = m.remove_row(0);
				m = m.remove_column(0);
				p = p - 1;
				q = q - 1;
			},
			GoToInitial => (),
			Halt(e) => break
		};
	}
	invariants
}

