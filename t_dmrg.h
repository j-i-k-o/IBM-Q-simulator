#pragma once

#include "itensor/all.h"
#include "itensor/mps/siteset.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <iostream>

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


	class myMPS{
		private:
			myvector<Index> a;  // Links
			myvector<Index> s;  // Sites
			myvector<ITensor> M; // 
			ITensor Psi; // \Psi^{c c+1}
			size_t L; //system size
			size_t pos; // [1:L-1] Psi^{pos pos+1}

		public:
			myMPS(const myvector<std::pair<double, double>>& coeffs, const myvector<Index> site_ind = myvector<Index>())
				: a(coeffs.size()-1), s(site_ind), M(coeffs.size()), L(coeffs.size()), pos(1){
				//generate MPS (coeff[1]|Up>+coeff[2]|Dn>)^L
				//|coeff[1]|**2+|coeff[2]|**2 must be 1
				// initialize Link
				for(auto& elem : this->a){
					elem = Index("LinkInd", 1, Link);
				}
				// initialize Site
				if(s.size() == 0){
					for(size_t i=0; i<coeffs.size(); i++){
						s.push_back(Index("SiteInd", 2, Site));
					}
				}

				assert(s.size() == coeffs.size());
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

				//set Mcur
				Psi = M[pos]*M[pos+1];
				//the system has already been normalized.
			}

			myMPS(const myMPS& obj) = default;
			myMPS(myMPS&&) = default;

			inline void normalize(){
				this->Psi /= norm(Psi);
			}

			void primeAll(){
				for(auto& elem : this->s){
					elem = prime(elem);
				}

				for(auto& elem : this->a){
					elem = prime(elem);
				}

				for(auto& elem : this->M){
					elem = prime(elem);
				}

				Psi = prime(Psi);
			}

			void decomposePsi(Args args = Args::global()){
				//\Psi_{a[i-1] a[i+1]}^{s[i], s[i+1]} -> M_{a[i-1], a[i]}^{s[i]}M_{a[i], a[i+1]}^{s[i+1]}
				ITensor U,S,V;
				if(this->pos == 1)
					U = ITensor(s[pos]);
				else
					U = ITensor(s[pos], a[pos-1]);

				Spectrum spec = svd(Psi, U, S, V, args);
					// M[i] = U, \Psi^{s[i+1] s[i+2]} = SV^{\dag}M[i+2]
				a[pos] = commonIndex(U,S,Link);
				S /= norm(S); //normalize
				this->M[pos] = U;
				this->M[pos+1] = S*V;
			}

			enum struct Direction{Right, Left};
			Spectrum shift_to(Direction dir, const Args& args = Args::global()){
				// shift \Psi to next site
				// returns current position (size_t)
				Spectrum spec;

				if(dir == Direction::Right){
					//sweep to right
					assert(pos < L-1);
					size_t i = pos;
					ITensor U,S,V;
					if(i == 1)
						U = ITensor(s[i]);
					else
						U = ITensor(s[i], a[i-1]);

					//M_{a[i-1] a[i+1]}^{s[i], s[i+1]} = U_{a[i-1] a[i]}^{s[i]}S_{a[i] a[i]}V_{a[i] a[i+1]}^{s[i+1]}
					spec = svd(Psi, U, S, V, args);
					// M[i] = U, \Psi^{s[i+1] s[i+2]} = SV^{\dag}M[i+2]
					a[i] = commonIndex(U,S,Link);
					M[i] = U;
					Psi = S*V*M[i+2];
					this->normalize();

					pos++;
					return spec;
				}
				else{
					//sweep to left
					assert(pos > 1);
					size_t i = pos;
					ITensor U,S,V;
					if(i == L-1)
						V = ITensor(s[i+1]);
					else
						V = ITensor(s[i+1], a[i+1]);

					//M_{a[i-1] a[i+1]}^{s[i], s[i+1]} = U_{a[i-1] a[i]}^{s[i]}S_{a[i] a[i]}V_{a[i] a[i+1]}^{s[i+1]}
					spec = svd(Psi, U, S, V, args);
					// M[i+1] = V, \Psi^{s[i-1] s[i]} = M[i-1]US
					a[i] = commonIndex(S,V,Link);
					M[i+1] = V;
					Psi = M[i-1]*U*S;
					this->normalize();

					pos--;
					return spec;
				}
			}

			inline size_t size() const{
				return L;
			}

			inline const ITensor& Mref_c(size_t i) const{
				assert(1<=i && i<=L);
				return M[i];
			}

			inline const myvector<ITensor>& Mref_c() const{
				return M;
			}

			inline const ITensor& Psi_c() const{
				return this->Psi;
			}

			inline const Index& site_c(size_t i) const{
				assert(1<=i && i<=L);
				return s[i];
			}

			inline const Index& link_c(size_t i) const{
				assert(1<=i && i<=L-1);
				return a[i];
			}

			inline const myvector<Index> site_c() const{
				return s;
			}

			inline const myvector<Index> link_c() const{
				return a;
			}
			inline const size_t get_pos() const{
				return this->pos;
			}

			void move_pos_to(size_t to, const Args& args = Args::global()){
				if(to < this->pos){
					//move left
					while(to != this->get_pos()){
						this->shift_to(myMPS::Direction::Left, args);
					}
				}
				if(to > this->pos){
					//move right
					while(to != this->get_pos()){
						this->shift_to(myMPS::Direction::Right, args);
					}
				}
			}

			void apply(const ITensor& Op){
				//apply operator (Op^{s[pos] s[pos+1]}) to wavefunction
				assert(Op.inds().size() == 4);
				for(auto&& elem : Op.inds()){
					//check if 
					assert(elem == s[pos] || elem == s[pos+1] || elem == prime(s[pos]) || elem == prime(s[pos+1]));
				}
				this->Psi = Op * prime(Psi,s[pos],s[pos+1]);
			}
	};

	//ITensor Hamiltonian
	inline const ITensor id(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(1), 1);
		ret.set(s(2), prime(s)(2), 1);
		return ret;
	}

	inline const ITensor sx(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(2), 1);
		ret.set(s(2), prime(s)(1), 1);
		return ret;
	}

	inline const ITensor sy(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(2), -1_i);
		ret.set(s(2), prime(s)(1), 1_i);
		return ret;
	}

	inline const ITensor sz(const Index& s){
		ITensor ret(s, prime(s));
		ret.set(s(1), prime(s)(1), 1);
		ret.set(s(2), prime(s)(2), -1);
		return ret;
	}
	inline const ITensor locH(double J1, double h1, double h2, double G1, double G2, double A, double B, const Index& s1, const Index& s2){
		return B * (J1*sz(s1)*sz(s2) + h1*sz(s1)*id(s2) + h2*id(s1)*sz(s2)) + A * (G1*sx(s1)*id(s2) + G2*id(s1)*sx(s2));
	}

	inline void printPsi(const myMPS& mps){
#ifdef DEBUG
			std::cout << "pos = " << mps.get_pos() << std::endl;
			PrintData(mps.Psi_c());
#endif
	}

	inline void printM(const myMPS& mps){
#ifdef DEBUG
		for(int i=1; i<=mps.size(); i++){
			std::cout << "i = " << i << std::endl;
			PrintData(mps.Mref_c(i));
		}
#endif
	}

	inline void printSpec(const Spectrum& spec){
#ifdef DEBUG
			PrintData(spec);
#endif
	}

	void sweep_Hodd(
			myMPS::Direction dir,
			myMPS& mps,
			const myvector<double>& J,
			const myvector<double>& h,
			const myvector<double>& G,
			const std::function<double(double)>& A, //X
			const std::function<double(double)>& B, //Z
			double time_arg,
			Cplx coef,
			const Args& args = Args::global()
			){
		Spectrum spec;

#ifdef DEBUG
		if(dir == myMPS::Direction::Right)
			std::cout << "======= sweep to right =======" << std::endl;
		else
			std::cout << "======= sweep to left =======" << std::endl;
#endif

		if(dir == myMPS::Direction::Right){
			//sweep to right
			//
			mps.move_pos_to(1, args);

			for(int pos=1; pos<=mps.size()-1; pos+=2){
				assert(mps.get_pos() == pos);
				mps.apply(
						expHermitian(locH(J[pos], h[pos], h[pos+1], G[pos], G[pos+1], A(time_arg), B(time_arg), mps.site_c(pos), mps.site_c(pos+1)), coef)
						);

				printM(mps);
				printPsi(mps);

				if(pos < mps.size()-1){
					//sweep twice if pos < mps.size()-1
					spec = mps.shift_to(myMPS::Direction::Right, args);
					printSpec(spec);
					spec = mps.shift_to(myMPS::Direction::Right, args);
					printSpec(spec);
				}
			}
		}
		else{
			//sweep to left
			//
			mps.move_pos_to(mps.size()-1, args);
			for(int pos=mps.size()-1; pos>=1; pos-=2){
				assert(mps.get_pos() == pos);
				mps.apply(
						expHermitian(locH(J[pos], h[pos], h[pos+1], G[pos], G[pos+1], A(time_arg), B(time_arg), mps.site_c(pos), mps.site_c(pos+1)), coef)
						);

				printM(mps);
				printPsi(mps);

				if(pos > 1){
					//sweep twice if pos > 1
					spec = mps.shift_to(myMPS::Direction::Left, args);
					printSpec(spec);
					spec = mps.shift_to(myMPS::Direction::Left, args);
					printSpec(spec);
				}
			}
		}
	}

	void sweep_Heven(
			myMPS::Direction dir,
			myMPS& mps,
			const myvector<double>& J,
			//const myvector<double>& h,
			//const myvector<double>& G,
			const std::function<double(double)>& A, //X
			const std::function<double(double)>& B, //Z
			double time_arg,
			Cplx coef,
			const Args& args = Args::global()
			){
		Spectrum spec;

#ifdef DEBUG
		if(dir == myMPS::Direction::Right)
			std::cout << "======= sweep to right =======" << std::endl;
		else
			std::cout << "======= sweep to left =======" << std::endl;
#endif

		if(dir == myMPS::Direction::Right){
			//sweep to right
			//
			mps.move_pos_to(2, args);

			for(int pos=2; pos<=mps.size()-2; pos+=2){
				assert(mps.get_pos() == pos);
				mps.apply(
						expHermitian(locH(J[pos], 0, 0, 0, 0, A(time_arg), B(time_arg), mps.site_c(pos), mps.site_c(pos+1)), coef)
						);

				printM(mps);
				printPsi(mps);

				if(pos < mps.size()-2){
					//sweep twice if pos < mps.size()-2
					spec = mps.shift_to(myMPS::Direction::Right, args);
					printSpec(spec);
					spec = mps.shift_to(myMPS::Direction::Right, args);
					printSpec(spec);
				}
			}
		}
		else{
			//sweep to left 
			//
			mps.move_pos_to(mps.size()-2, args);
			for(int pos=mps.size()-2; pos>=2; pos-=2){
				assert(mps.get_pos() == pos);
				mps.apply(
						expHermitian(locH(J[pos], 0, 0, 0, 0, A(time_arg), B(time_arg), mps.site_c(pos), mps.site_c(pos+1)), coef)
						);

				printM(mps);
				printPsi(mps);

				if(pos > 2){
					//sweep twice if pos > 2
					spec = mps.shift_to(myMPS::Direction::Left, args);
					printSpec(spec);
					spec = mps.shift_to(myMPS::Direction::Left, args);
					printSpec(spec);
				}
			}
		}
	}

	void timeevolve(
			myMPS& mps,
			const myvector<double>& J,
			const myvector<double>& h,
			const myvector<double>& G,
			const std::function<double(double)>& A, //X
			const std::function<double(double)>& B, //Z
			double t,
			double dt,
			const Args& args = Args::global()
			){
		assert(mps.size()%2 == 0);
		assert(J.size() == mps.size()-1);
		assert(h.size() == mps.size());
		assert(G.size() == mps.size());

		double s2 = 1/(4-pow(4,1.0/3.0));

		////for more details, see Hatano and Suzuki (2005)
		//
		//// e ^ -i/2 s2 dt Hodd (t+     s2/2 dt)
		//sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt*s2/2., -s2 * dt * 1_i/2., args);
		//// e ^ -i   s2 dt Heven(t+     s2/2 dt)                                     
		//sweep_Heven(myMPS::Direction::Left, mps, J,       A, B, t+dt*s2/2., -s2 * dt * 1_i   , args);
		//// e ^ -i/2 s2 dt Hodd (t+     s2/2 dt)                                     
		//sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt*s2/2., -s2 * dt * 1_i/2., args);

		//mps.normalize();
		////
		//// e ^ -i/2 s2 dt Hodd (t+   3*s2/2 dt)
		//sweep_Hodd(myMPS::Direction::Left, mps, J, h, G, A, B, t+3*dt*s2/2., -s2 * dt * 1_i/2., args);
		//// e ^ -i   s2 dt Heven(t+   3*s2/2 dt)                                      
		//sweep_Heven(myMPS::Direction::Right, mps, J,     A, B, t+3*dt*s2/2., -s2 * dt * 1_i   , args);
		//// e ^ -i/2 s2 dt Hodd (t+   3*s2/2 dt)                                      
		//sweep_Hodd(myMPS::Direction::Left, mps, J, h, G, A, B, t+3*dt*s2/2., -s2 * dt * 1_i/2., args);

		//mps.normalize();
		////  
		//// e ^ -i/2 s2 dt Hodd (t+      1/2 dt)
		//sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt/2., -s2 * dt * 1_i/2., args);
		//// e ^ -i   s2 dt Heven(t+      1/2 dt)                                  
		//sweep_Heven(myMPS::Direction::Left, mps, J,       A, B, t+dt/2., -s2 * dt * 1_i   , args);
		//// e ^ -i/2 s2 dt Hodd (t+      1/2 dt)                                  
		//sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt/2., -s2 * dt * 1_i/2., args);

		//mps.normalize();
		////
		//// e ^ -i/2 s2 dt Hodd (t+dt-3*s2/2 dt)
		//sweep_Hodd(myMPS::Direction::Left, mps, J, h, G, A, B, t+dt-3*dt*s2/2., -s2 * dt * 1_i/2., args);
		//// e ^ -i   s2 dt Heven(t+dt-3*s2/2 dt)                                         
		//sweep_Heven(myMPS::Direction::Right, mps, J,     A, B, t+dt-3*dt*s2/2., -s2 * dt * 1_i   , args);
		//// e ^ -i/2 s2 dt Hodd (t+dt-3*s2/2 dt)                                         
		//sweep_Hodd(myMPS::Direction::Left, mps, J, h, G, A, B, t+dt-3*dt*s2/2., -s2 * dt * 1_i/2., args);

		//mps.normalize();
		////  
		//// e ^ -i/2 s2 dt Hodd (t+dt-  s2/2 dt)
		//sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt-dt*s2/2., -s2 * dt * 1_i/2., args);
		//// e ^ -i   s2 dt Heven(t+dt-  s2/2 dt)                                        
		//sweep_Heven(myMPS::Direction::Left, mps, J,       A, B, t+dt-dt*s2/2., -s2 * dt * 1_i   , args);
		//// e ^ -i/2 s2 dt Hodd (t+dt-  s2/2 dt)                                        
		//sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt-dt*s2/2., -s2 * dt * 1_i/2., args);

		//mps.normalize();

		//// e ^ -i dt Hodd (t+     dt)
		//sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt, -dt * 1_i, args);
		//// e ^ -i dt Heven(t+     dt)                                     
		//sweep_Heven(myMPS::Direction::Left, mps, J,       A, B, t+dt, -dt * 1_i, args);
		// return
		// e ^ -i/2 dt Hodd (t+     dt/2)
		sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt/2., -dt * 1_i/2., args);
		// e ^ -i dt Heven(t+     dt)                                     
		sweep_Heven(myMPS::Direction::Left, mps, J,       A, B, t+dt/2., -dt * 1_i, args);
		// e ^ -i/2 dt Hodd (t+     dt/2)
		sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, t+dt/2., -dt * 1_i/2., args);
	}

	Cplx overlap(const myMPS& t_mps1, const myvector<ITensor>& op, const myMPS& t_mps2, const Args& args = Args::global()){
		assert(op.size() == t_mps1.size());
		assert(op.size() == t_mps2.size());

		//create copy
		myMPS mps1 = t_mps1;
		myMPS mps2 = t_mps2;
		//move to middle

		mps1.move_pos_to(mps1.size()/2, args);
		mps2.move_pos_to(mps2.size()/2, args);
		mps1.move_pos_to(mps1.size()/2, args);
		mps2.move_pos_to(mps2.size()/2, args);
		mps1.decomposePsi(args);
		mps2.decomposePsi(args);

		printM(mps1);
		printPsi(mps1);
		printM(mps2);
		printPsi(mps2);

		mps2.primeAll();

		ITensor ret_t = dag(mps1.Mref_c(1))*op[1]*mps2.Mref_c(1);
		for(size_t i=2; i<= mps1.size(); i++){
			//reduction
			ret_t = dag(mps1.Mref_c(i))*op[i]*ret_t*mps2.Mref_c(i);
		}
		
		return ret_t.cplx();
	}




	

};
