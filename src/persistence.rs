use crate::topology::*;
use crate::topology::homology::Homology;
use kiddo::KdTree;
use kiddo::Manhattan;
use core::array::from_fn;


pub struct VietorisRips<const D: usize, const T: usize>{
    ε: f64,
    pts: usize,
    tree: KdTree<f64,D>,
    skeleton: Vec<Vec<bool>>,
    raw_complex: Vec<Simplex>,
}
impl<const D: usize, const T: usize> VietorisRips<D,T>{
    pub fn init(cloud: Vec<[f64;D]>) -> VietorisRips<{D},{T}>{
        assert!(D <= 3 && D >= 2,"supported dimensions: 2, 3.");
        let mut i = 0;
        let pts = cloud.len();

        let mut tree: KdTree<f64, {D}> = KdTree::new();
        let skeleton: Vec<Vec<bool>> = vec![vec![false;pts];pts];
        let mut raw_complex: Vec<Simplex> = Vec::new();

        for array in cloud.iter(){
            tree.add(&array,(pts - i - 1) as u64);
            raw_complex.push(Simplex{
                value: vec![i as usize],
                dimension: 0
            });
            i += 1;
        }
        VietorisRips{
            ε: 0.0,
            pts,
            tree,
            skeleton,
            raw_complex
        }
    }
    pub fn add_2_simplex(&mut self, triangles: &mut Vec<[usize;3]>) -> bool{
        let mut added_2_simplex = false;
        for a0 in 0..self.pts{
            for a1 in a0..self.pts{
                if self.skeleton[a0][a1]{
                    for a2 in a1..self.pts{
                        if self.skeleton[a1][a2] && self.skeleton[a0][a2] && !triangles.contains(&[a0,a1,a2]){
                            triangles.push([a0,a1,a2]);
                            self.raw_complex.push(Simplex{
                                value: vec![a0,a1,a2],
                                dimension: 2
                            });
                            added_2_simplex = true;
                        }
                    }
                }
            }
        }
        added_2_simplex
    }
    pub fn add_3_simplex(&mut self, tetrahedrons: &mut Vec<[usize;4]>){
        for a0 in 0..self.pts{
            for a1 in a0..self.pts{
                if self.skeleton[a0][a1]{
                    for a2 in a1..self.pts{
                        if self.skeleton[a1][a2]{
                            for a3 in a2..self.pts{
                                if self.skeleton[a0][a2] && self.skeleton[a0][a3] && self.skeleton[a1][a3] && self.skeleton[a2][a3] && !tetrahedrons.contains(&[a0,a1,a2,a3]){
                                    tetrahedrons.push([a0,a1,a2,a3]);
                                    self.raw_complex.push(Simplex{
                                        value: vec![a0,a1,a2,a3],
                                        dimension: 3
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    pub fn analyze(&mut self, ring: Ring, end: f64, step: f64) -> Vec<Homology<{T}>>{

        let mut first_iteration = true;
        let mut added_1_simplex = false;
        let mut triangles: Vec<[usize;3]> = Vec::new();
        let mut tetrahedrons: Vec<[usize;4]> = Vec::new();

        let mut current_homology = Homology{
            start: 0.0,
            end: 0.0,
            ranks: [0;T],
            torsions: from_fn(|_| vec![])
        };
        let mut data: Vec<Homology<T>> = Vec::new();

        self.ε = 0.0;
        while self.ε <= end{
            println!("ε={}", self.ε);
            for (i,pos) in self.tree.iter(){
                let nbhd = self.tree.within::<Manhattan>(&pos,self.ε);
                for point in nbhd{
                    if point.item > i && point.distance >= self.ε - step{
                        self.raw_complex.push(Simplex{
                            value: vec![i as usize,point.item as usize],
                            dimension: 1
                        });

                        self.skeleton[i as usize][point.item as usize] = true;
                        self.skeleton[point.item as usize][i as usize] = true;

                        added_1_simplex = true;
                    }
                }
            }

            if added_1_simplex{ 
                if self.add_2_simplex(&mut triangles) && D == 3{
                    self.add_3_simplex(&mut tetrahedrons);
                }
            }
            if added_1_simplex || first_iteration{
                self.update_homology(SimplicialComplex::from(self.raw_complex.clone(), ring.clone()), &mut current_homology, &mut data, first_iteration);
            }

            first_iteration = false;
            added_1_simplex = false;
            self.ε += step;
        }
        data.push(Homology{
            start: current_homology.start,
            end: self.ε,
            ranks: current_homology.ranks,
            torsions: current_homology.torsions.clone()
        });
        data
    }

    pub fn update_homology(&mut self, vr: SimplicialComplex, current_homology: &mut Homology<T>, data: &mut Vec<Homology<T>>, first_iteration: bool){
        let ranks_and_torsions: [(usize,Torsion);T] = from_fn(|k| vr.H(k));
        let new_ranks = ranks_and_torsions.clone().map(|c| c.0);
        let new_torsions = ranks_and_torsions.map(|c| c.1);
        if first_iteration{
            current_homology.start = self.ε;
            current_homology.ranks = new_ranks;
            current_homology.torsions = new_torsions.clone();
        }
        if current_homology.ranks != new_ranks || current_homology.torsions != new_torsions{

            data.push(Homology{
                start: current_homology.start,
                end: self.ε,
                ranks: current_homology.ranks,
                torsions: current_homology.torsions.clone()
            });
            current_homology.start = self.ε;
            current_homology.ranks = new_ranks;
            current_homology.torsions = new_torsions;
        }
    }
}
