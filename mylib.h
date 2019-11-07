#pragma once

#include "itensor/all.h"
#include "itensor/mps/siteset.h"
#include <vector>
#include <utility>
#include <type_traits>

namespace myLib{

	using namespace itensor;

	/******
	 * 1-indexed vector
	 ******/

	template<typename _Tp, typename _Alloc = std::allocator<_Tp>>
		class myvector : public std::vector<_Tp, _Alloc>
	{
		using _Base = std::vector<_Tp, _Alloc>;

		public:
			myvector() : _Base(){}

			explicit
				myvector(size_t __n, const _Alloc& __a = _Alloc()) : _Base(__n, __a){}

			myvector(const myvector& __x) : _Base(__x){}
			myvector(myvector&& __x) noexcept : _Base(__x){}

			myvector(const _Base& __x) : _Base(__x){}
			myvector(_Base&& __x) noexcept : _Base(__x){}

			myvector& operator=(const myvector& __x){
				this->clear();
				for(auto& elem : __x){
					this->push_back(elem);
				}
				return *this;
			}

			myvector& operator=(const _Base& __x){
				return this->operator=(__x);
			}

			//index from one
			
			_Tp& operator[](size_t __n){
				return _Base::operator[](__n - 1);
			}

			const _Tp& operator[](size_t __n) const{
				return _Base::operator[](__n - 1);
			}
	};


	// MPS class
	
	class myMPS;
	class myMPO;
	class LadderMPO;
	
	class myMPS{
		private:
			myvector<Index> a;  // Links
			myvector<Index> s;  // Sites
			myvector<ITensor> M; // ITensors
			size_t L; //system size

				
			size_t ortho_pos; // 0: not orthogonalized, >0: where the center is.

			void init_MPS(size_t L, long link_size, bool is_initSite){
				//generate Link Indices
				for(Index& elem : a){
					elem = Index("LinkInd", link_size, Link);
				}
				//generate Site Indices
				if(is_initSite){
					for(Index& elem : s){
						elem = Index("SiteInd", 2, Site);
					}
				}
				//generate ITensors
				for(auto l : range1(L)){
					if(l == 1){
						M[l] = randomTensor(s[l], a[l]);
					}
					else if(l == L){
						M[l] = randomTensor(s[l], a[l-1]);
					}
					else{
						M[l] = randomTensor(s[l], a[l-1], a[l]);
					}
				}
			}


			void init_MPS(const myvector<ITensor>& refM){
				//generate myMPS class from refM
				// refM[i] -> ITensor(Site, Link, (Link))
				this->M = refM;
				//set indices
				set_indices();
			}

			void init_MPS(const myvector<std::pair<double, double>>& coeffs){
				//generate MPS (coeff[1st]|Up>+coeff[2nd]|Dn>)^L
				// initialize Link
				for(auto& elem : this->a){
					elem = Index("LinkInd", 1, Link);
				}
				// initialize Site
				for(auto& elem : this->s){
					elem = Index("SiteInd", 2, Site);
				}
				// initialize Tensors
				for(auto l : range1(L)){
					if(l == 1){
						M[l] = ITensor(s[l], a[l]);
						M[l].set(a[l](1), s[l](1), coeffs[l].first);
						M[l].set(a[l](1), s[l](2), coeffs[l].second);
					}
					else if(l == L){
						M[l] = ITensor(s[l], a[l-1]);
						M[l].set(a[l-1](1), s[l](1), coeffs[l].first);
						M[l].set(a[l-1](1), s[l](2), coeffs[l].second);
					}
					else{
						M[l] = ITensor(s[l], a[l-1], a[l]);
						M[l].set(a[l-1](1), a[l](1), s[l](1), coeffs[l].first);
						M[l].set(a[l-1](1), a[l](1), s[l](2), coeffs[l].second);
						
					}
				}
			}

		public:

			//generate random MPS
			myMPS(size_t L, long link_size = 1)
				: a(L-1), s(L), M(L), L(L), ortho_pos(0){
					// sites will be initialized
					init_MPS(L, link_size, true);
			}

			myMPS(const myvector<Index>& sites, long link_size = 1)
				: a(sites.size()-1), s(sites), M(sites.size()), L(sites.size()), ortho_pos(0){
					// sites will *not* be initialized
					init_MPS(L, link_size, false);
			}

			myMPS(const myvector<ITensor>& tensors)
				: a(tensors.size()-1), s(tensors.size()), M(tensors.size()), L(tensors.size()), ortho_pos(0){
					init_MPS(tensors);
		 	}

			myMPS(const myvector<std::pair<double, double>>& coeffs)
				: a(coeffs.size()-1), s(coeffs.size()), M(coeffs.size()), L(coeffs.size()), ortho_pos(0){
				//generate MPS (coeff[1]|Up>+coeff[2]|Dn>)^L
				init_MPS(coeffs);
			}

			void set_indices(){
				// set a[i] and s[i] from M
				//set indices
				for(size_t i : range1(M.size())){
					//set Site Index
					for(auto& elem : M[i].inds()){
						if(elem.type() == Site){
							this->s[i] = elem;
							break;
						}
					}
					//set Link Index
					if(i != L)
						this->a[i] = commonIndex(M[i], M[i+1]);
				}

				for(auto& elem : this->s){
					assert(elem.type() == Site);
				}
				for(auto& elem : this->a){
					assert(elem.type() == Link);
				}
			}

			Spectrum mix_canonical(size_t i, const Args& args = Args::global()){
				// M[1:i-1] -> left-normalized matrices (often referred to as A)
				// M[i+1:L] -> right-normalized matrices (often referred to as B)
				// M[i]     -> unnormalized matrices (often referred to as \Psi)
				// |\psi> = AAAAA...A\PsiB...BBBBB
				Spectrum spec;
				assert(0 <= ortho_pos && ortho_pos <= L);
				if(ortho_pos == 0){
					// not have been orthogonalized
					ortho_pos = 1;
					// sweep
					while(ortho_pos != L)
						shift_to(Direction::Right, args);
					while(ortho_pos != i)
						spec = shift_to(Direction::Left, args);
				}
				else{
					// already orthogonalized
					if(i > ortho_pos){
						while(ortho_pos != i)
							shift_to(Direction::Right, args);
					}
					if(i < ortho_pos){
						while(ortho_pos != i)
							spec = shift_to(Direction::Right, args);
					}
				}
				return spec;
			}

			enum struct Direction{Right, Left};
			Spectrum shift_to(Direction dir, const Args& args = Args::global()){
				// shift \Psi to next site
				// returns current position (size_t)
				// assume that MPS has been orthogonalized
				assert(ortho_pos != 0);
				Spectrum spec;

				if(dir == Direction::Right){
					//sweep to right
					if(ortho_pos == L){
						return spec;
					}
					size_t i = ortho_pos;
					ITensor Psi = M[i];
					ITensor U,S,V;
					if(i == 1)
						U = ITensor(s[i]);
					else
						U = ITensor(s[i], a[i-1]);

					spec = svd(Psi, U, S, V, args);
					// A[i] = U, \Psi[i+1] = SV^{\dag}M[i+1]
					a[i] = commonIndex(U,S,Link);
					M[i] = U;
					M[i+1] = S*V*M[i+1];

					ortho_pos++;
					return spec;

				}
				else{
					//sweep to left
					if(ortho_pos == 1){
						return spec;
					}
					size_t i = ortho_pos;
					ITensor Psi = M[i];
					ITensor U,S,V;
					if(i == L)
						V = ITensor(s[i]);
					else
						V = ITensor(s[i], a[i]);

					spec = svd(Psi, U, S, V, args);
					// B[i] = V, \Psi[i-1] = M[i-1]US
					a[i-1] = commonIndex(S,V,Link);
					M[i] = V;
					M[i-1] = M[i-1]*U*S;

					ortho_pos--;
					return spec;
				}
			}

			
			Real norm() const{
				assert(ortho_pos != 0);
				//M[i](a[i-1], s[i], a[i]) -> U(a[i-1], s[i], k)S(k, k)V(k, a[i])
				size_t i = ortho_pos;
				return itensor::norm(M[i]); 
			}
			

			
			Real normalize(){
				assert(ortho_pos != 0);
				//M[i](a[i-1], s[i], a[i]) -> U(a[i-1], s[i], k)S(k, k)V(k, a[i])
				size_t i = ortho_pos;
				Real nrm = itensor::norm(M[i]); 
				M[i] *= 1./nrm;
				return nrm;
			}
			

			inline size_t position() const{
				return ortho_pos;
			}


			//Tensors
			const myvector<ITensor>& Mref_c() const{
				return this->M;
			}

			const ITensor& Mref_c(size_t i) const{
				assert(1<=i && i<=L);
				return M[i];
			}
			
			ITensor& Mref(size_t i){
				assert(1<=i && i<=L);
				// MPS is not guaranteed to be orthogonalized.
				// ortho_pos is set to zero.
				ortho_pos = 0;
				return M[i];
			}

			//indices
			inline const myvector<Index>& links_ref() const{
				return this->a;
			}

			inline const myvector<Index>& sites_ref() const{
				return this->s;
			}

			//set prime
			void mapprime_Site(int newlevel = 1){
				for(size_t i : range1(L)){
					M[i].mapprime(s[i].primeLevel(), newlevel, Site);
					s[i].primeLevel(newlevel);
				}
			}

			inline void mapnoprime_Site(){
				mapprime_Site(0);
			}

	};


	//fundamental operators
	// pauli_x
	inline Real s_x(size_t s1, size_t s2){
		if(s1 == 1 && s2 == 2) return 1;
		else if(s1 == 2 && s2 == 1) return 1;
		else return 0;
	}
	// pauli_y
	inline Cplx s_y(size_t s1, size_t s2){
		if(s1 == 1 && s2 == 2) return -1_i;
		else if(s1 == 2 && s2 == 1) return 1_i;
		else return 0;
	}
	// pauli_z
	inline Real s_z(size_t s1, size_t s2){
		if(s1 == 1 && s2 == 1) return 1;
		else if(s1 == 2 && s2 == 2) return -1;
		else return 0;
	}
	// id 
	inline Real id(size_t s1, size_t s2){
		if(s1 == 1 && s2 == 1) return 1;
		else if(s1 == 2 && s2 == 2) return 1;
		else return 0;
	}

	//MPO
	class myMPO{
		protected:
			myvector<Index> a; // Links
			myvector<Index> s; // Sites
			myvector<ITensor> W; // ITensors
			size_t L; //system size

			void init_MPO(const myvector<ITensor>& refW){
				//initialize general MPO with vector W
				this->W = refW;
				set_indices();
			}

		public:

			myMPO(const myvector<Index>& a, const myvector<Index>& s, const myvector<ITensor>& W, size_t L) : a(a), s(s), W(W), L(L){}

			myMPO(const myvector<ITensor>& refW) : a(refW.size()-1), s(refW.size()), L(refW.size()){
				init_MPO(refW);
			}

			void set_indices(){
				//set indices from W
				for(size_t i : range1(W.size())){
					//set Site Index
					for(auto& elem : W[i].inds()){
						if(elem.type() == Site){
							this->s[i] = elem;
							break;
						}
					}
					//set Link Index
					if(i != L)
						this->a[i] = commonIndex(W[i], W[i+1]);
				}

				for(auto& elem : this->s){
					assert(elem.type() == Site);
				}
				for(auto& elem : this->a){
					assert(elem.type() == Link);
				}
			}

			//indices
			inline const myvector<Index>& links_ref() const{
				return this->a;
			}

			inline const myvector<Index>& sites_ref() const{
				return this->s;
			}

			//Tensors
			const myvector<ITensor>& Wref_c() const{
				return this->W;
			}
			
			const ITensor& Wref_c(size_t i) const{
				assert(1<=i && i<=L);
				return W[i];
			}
			
			ITensor& Wref(size_t i){
				assert(1<=i && i<=L);
				return W[i];
			}

			void prime_Site(int inc = 1){
				for(size_t i : range1(L)){
					W[i].prime(Site, inc);
					s[i].prime(inc);
				}
			}

	};


	// MPO for quasi-1d ladder Hamiltonian
	class LadderMPO : public myMPO{
		private:
			//J1 -> three body interaction
			//J2 -> two body NNN
			//J1.size() == systemsize - 2
			void init_LadderMPO(const myvector<double>& J1, const myvector<double>& J2, double Gamma, bool is_initSite){
				size_t L = this->L; //system size
				assert(J1.size() == J2.size());
				//generate Link Indices
				for(Index& elem : a){
					elem = Index("LinkInd", 4, Link);
				}
				//generate Site Indices
				if(is_initSite){
					for(Index& elem : s){
						elem = Index("SiteInd", 2, Site);
					}
				}
				//generate ITensors
				//
				// 1 < l < L
				//
				// W=
				//   [  1   0          0   0  ]
				//   [  sz  0          0   0  ]
				//   [  1   (J1sz+J2)  0   0  ]
				//   [  sx  0          sz  1  ]
				//
				// l = 1
				//
				// W=
				//   [  sx  0          sz  1  ]
				//
				// l = L
				//
				// W=
				//   [  1   ]
				//   [  sz  ]
				//   [  1   ]
				//   [  sx  ]
				//


				//TODO: more efficient (and sophisticated) method?
				for(auto l : range1(L)){
					if(l == 1){
						W[l] = ITensor(s[l], prime(s[l]), a[l]);
						for(auto i : range1(2))
						for(auto j : range1(2))
						{
								W[l].set(a[l](1), s[l](i), prime(s[l])(j), Gamma*s_x(i, j));
								W[l].set(a[l](2), s[l](i), prime(s[l])(j), 0);
								W[l].set(a[l](3), s[l](i), prime(s[l])(j), s_z(i, j));
								W[l].set(a[l](4), s[l](i), prime(s[l])(j), id(i, j));
						}	
					}
					else if(l == L){
						W[l] = ITensor(s[l], prime(s[l]), a[l-1]);
						for(auto i : range1(2))
						for(auto j : range1(2))
						{
								W[l].set(a[l-1](1), s[l](i), prime(s[l])(j), id(i, j));
								W[l].set(a[l-1](2), s[l](i), prime(s[l])(j), s_z(i, j));
								W[l].set(a[l-1](3), s[l](i), prime(s[l])(j), 0);
								W[l].set(a[l-1](4), s[l](i), prime(s[l])(j), Gamma*s_x(i, j));
						}	
					}
					else{
						W[l] = ITensor(s[l], prime(s[l]), a[l-1], a[l]);
						for(auto i : range1(2))
						for(auto j : range1(2))
						{
								W[l].set(a[l-1](1), a[l](1), s[l](i), prime(s[l])(j), id(i, j));
								W[l].set(a[l-1](2), a[l](1), s[l](i), prime(s[l])(j), s_z(i, j));
								W[l].set(a[l-1](3), a[l](2), s[l](i), prime(s[l])(j), J1[l-1]*s_z(i, j)+J2[l-1]*id(i, j));
								W[l].set(a[l-1](4), a[l](1), s[l](i), prime(s[l])(j), Gamma*s_x(i, j));
								W[l].set(a[l-1](4), a[l](3), s[l](i), prime(s[l])(j), s_z(i, j));
								W[l].set(a[l-1](4), a[l](4), s[l](i), prime(s[l])(j), id(i, j));
						}
					}
				}
			}
		public:
			LadderMPO(const myvector<double>& J1, const myvector<double>& J2, double Gamma) 
				: myMPO(myvector<Index>(J1.size()+1), myvector<Index>(J1.size()+2), myvector<ITensor>(J1.size()+2), J1.size()+2){
					// sites will be initialized
				init_LadderMPO(J1, J2, Gamma, true);
			}

			//with explicitly determined sites
			LadderMPO(const myvector<Index>& sites, const myvector<double>& J1, const myvector<double>& J2, double Gamma)
				: myMPO(myvector<Index>(J1.size()+1), sites, myvector<ITensor>(J1.size()+2), J1.size()+2){
				assert(sites.size() == J1.size() + 2);
				assert(sites[1].m() == 2);
					// sites will *not* be initialized
				init_LadderMPO(J1, J2, Gamma, false);
			}

	};


	/**********************
	 * Utility functions
	 **********************/

	//TODO: define operator= operator=(&&)

	//overlap <l|r>
	Real overlap(const myMPS& mps_l, const myMPS& mps_r);

	//add
	myMPS add(const myMPS& mps_l, const myMPS& mps_r);

	//expectation value
	Real expect(const myMPO& mpo, const myMPS& mps);

	//expectation value(with local operator)
	Real expect(const ITensor& site_op, const myMPS& mps);

	//apply
	myMPS apply(const myMPO& mpo, const myMPS& mps);

	/***************************
	 * Optimization functions
	 ***************************/

	//compress
	//TODO: consider how to determine the tolerance
	myMPS compress(const myMPS& mps, const Args& args = Args::global());

	//minimize (contains DMRG)
	//TODO: consider how to determine the tolerance
	myMPS minimize(const myMPO& mpo, const myMPS& mps, const Args& args = Args::global());

	//TODO: time evolution



}
