#include "gtest/gtest.h"
#include "t_dmrg.h" 
using namespace itensor;
using namespace myLib;

class MyTest : public ::testing::Test{};

TEST_F(MyTest, sampletest){
	EXPECT_EQ(true, true);
}

TEST_F(MyTest, normcheck){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	myLib::myMPS mps(coeffs);
	auto Psi = mps.Psi_c();
	EXPECT_NEAR(norm(Psi), 1.0, 1e-8);

	Spectrum spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	auto Psi2 = mps.Psi_c();
	EXPECT_NEAR(norm(Psi2), 1.0, 1e-8);
}

TEST_F(MyTest, copyconstructorcheck){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}
	myLib::myMPS mps(coeffs);
	myLib::myMPS mps2 = mps;
	EXPECT_EQ(mps.site_c(), mps2.site_c());
	EXPECT_EQ(mps.link_c(), mps2.link_c());

	EXPECT_TRUE(hasindex(mps.Psi_c(), mps2.site_c(1)));
	EXPECT_TRUE(hasindex(mps.Psi_c(), mps2.site_c(2)));
	EXPECT_TRUE(hasindex(mps.Psi_c(), mps2.link_c(2)));
	EXPECT_EQ(mps.Psi_c().inds().size(), 3);
	EXPECT_EQ(mps2.Psi_c().inds().size(), 3);
}

TEST_F(MyTest, indexcheck){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	myLib::myMPS mps(coeffs); //1
	Spectrum spec = mps.shift_to(myMPS::Direction::Right); //2
	spec = mps.shift_to(myMPS::Direction::Right); //3
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right); //7
	spec = mps.shift_to(myMPS::Direction::Right); 
	spec = mps.shift_to(myMPS::Direction::Right); //9
	spec = mps.shift_to(myMPS::Direction::Left); 
	spec = mps.shift_to(myMPS::Direction::Left); 
	spec = mps.shift_to(myMPS::Direction::Left); 
	spec = mps.shift_to(myMPS::Direction::Left); //5

	myvector<ITensor> M = mps.Mref_c();
	ITensor Psi = mps.Psi_c();
	myvector<Index> s = mps.site_c();
	myvector<Index> a = mps.link_c();

	IndexSet inds = Psi.inds();
	EXPECT_TRUE(hasindex(inds, s[5]));
	EXPECT_TRUE(hasindex(inds, s[6]));
	EXPECT_TRUE(hasindex(inds, a[4]));
	EXPECT_TRUE(hasindex(inds, a[6]));
	EXPECT_EQ(inds.size(), 4);

	EXPECT_TRUE(hasindex(M[1].inds(), s[1]));
	EXPECT_TRUE(hasindex(M[1].inds(), a[1]));
	EXPECT_EQ(M[1].inds().size(), 2);

	for(auto&& i : {2,3,4,7,8,9}){
		EXPECT_TRUE(hasindex(M[i].inds(), s[i]));
		EXPECT_TRUE(hasindex(M[i].inds(), a[i-1]));
		EXPECT_TRUE(hasindex(M[i].inds(), a[i]));
		EXPECT_EQ(M[i].inds().size(), 3);
	}
	
	EXPECT_TRUE(hasindex(M[10].inds(), s[10]));
	EXPECT_TRUE(hasindex(M[10].inds(), a[9]));
	EXPECT_EQ(M[10].inds().size(), 2);
}

TEST_F(MyTest, siteimmutablecheck){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	myLib::myMPS mps(coeffs);
	myvector<Index> s = mps.site_c();

	//sweep
	Spectrum spec;
	for(size_t i=1; i<9; i++){
		spec = mps.shift_to(myMPS::Direction::Right);
	}
	for(size_t i=1; i<9; i++){
		spec = mps.shift_to(myMPS::Direction::Left); 
	}
	for(size_t i=1; i<9; i++){
		spec = mps.shift_to(myMPS::Direction::Right);
	}
	for(size_t i=1; i<9; i++){
		spec = mps.shift_to(myMPS::Direction::Left); 
	}

	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);

	EXPECT_TRUE(hasindex(mps.Psi_c(), s[3]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[4]));

	for(auto&& i : {1,2,5,6,7,8,9,10}){
		EXPECT_TRUE(hasindex(mps.Mref_c(i), s[i]));
	}
}

TEST_F(MyTest, applycheck){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	myLib::myMPS mps(coeffs);

	Spectrum spec = mps.shift_to(myMPS::Direction::Right); //2
	spec = mps.shift_to(myMPS::Direction::Right); //3
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right); //6

	size_t i = mps.get_pos(); //6
	auto Op = ITensor(mps.site_c(i), mps.site_c(i+1), prime(mps.site_c(i)), prime(mps.site_c(i+1)));
	Op.set(mps.site_c(i)(1), mps.site_c(i+1)(1), prime(mps.site_c(i))(1), prime(mps.site_c(i+1))(1), 1);

	mps.apply(Op);

	myvector<Index> s = mps.site_c();
	myvector<Index> a = mps.link_c();

	EXPECT_TRUE(hasindex(mps.Psi_c(), s[i]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[i+1]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[i-1]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[i+1]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 4);
}

TEST_F(MyTest, antipatterncheck){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	myLib::myMPS mps(coeffs);
	myvector<Index> s = mps.site_c();
	myvector<Index> a = mps.link_c(); //this vector will be invalid!! do not do this

	Spectrum spec;
	spec = mps.shift_to(myMPS::Direction::Right); //2
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right);
	spec = mps.shift_to(myMPS::Direction::Right); 
	spec = mps.shift_to(myMPS::Direction::Right); //6
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[6]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[7]));
	EXPECT_FALSE(hasindex(mps.Psi_c(), a[5])); // this Index has been changed
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[7]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 4);
	EXPECT_EQ(mps.get_pos(), 6);

	a = mps.link_c(); //fetch link Index after shift_to operation completed
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[6]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[7]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[5]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[7]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 4);
	EXPECT_EQ(mps.get_pos(), 6);
}

TEST_F(MyTest, movepostocheck){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	myLib::myMPS mps(coeffs);
	//move
	myvector<Index> s = mps.site_c();
	mps.move_pos_to(7);
	myvector<Index> a = mps.link_c();
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[7]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[8]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[6]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[8]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 4);
	EXPECT_EQ(mps.get_pos(), 7);
	//move
	mps.move_pos_to(3);
	a = mps.link_c();
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[3]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[4]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[2]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[4]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 4);
	EXPECT_EQ(mps.get_pos(), 3);
	mps.move_pos_to(9);
	a = mps.link_c();
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[9]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[10]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[8]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 3);
	EXPECT_EQ(mps.get_pos(), 9);
	//move
	mps.move_pos_to(1);
	a = mps.link_c();
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[1]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[2]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[2]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 3);
	EXPECT_EQ(mps.get_pos(), 1);
	mps.move_pos_to(5);
	a = mps.link_c();
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[5]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[6]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[4]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[6]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 4);
	EXPECT_EQ(mps.get_pos(), 5);
}

TEST_F(MyTest, tebdtest){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	//initial state => ground state of -sx
	myLib::myMPS mps(coeffs);

	//initialize each coefficient
	myvector<double> J = std::vector<double>(9, -1);
	myvector<double> h = std::vector<double>(10, -1);
	myvector<double> G = std::vector<double>(10, -1);

	auto A = [](double t){return 1-t;};
	auto B = [](double t){return t;};

	double s2 = 1/(4-pow(4,1/3.));

	sweep_Hodd(myMPS::Direction::Right, mps, J, h, G, A, B, 0+0.01*s2/2., -s2 * 1_i/2., {"Cutoff", 1E-4});

	myvector<Index> s = mps.site_c();
	myvector<Index> a = mps.link_c();

	EXPECT_TRUE(hasindex(mps.Psi_c(), s[9]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[10]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[8]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 3);
	EXPECT_EQ(mps.get_pos(), 9);

	sweep_Heven(myMPS::Direction::Left, mps, J, A, B, 0+0.01*s2/2., -s2 * 1_i/2., {"Cutoff", 1E-4});
	a = mps.link_c();
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[2]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[3]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[1]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[3]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 4);
	EXPECT_EQ(mps.get_pos(), 2);

	sweep_Hodd(myMPS::Direction::Left, mps, J, h, G, A, B, 0+0.01*s2/2., -s2 * 1_i/2., {"Cutoff", 1E-4});
	a = mps.link_c();
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[1]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[2]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[2]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 3);
	EXPECT_EQ(mps.get_pos(), 1);

	sweep_Heven(myMPS::Direction::Right, mps, J, A, B, 0+0.01*s2/2., -s2 * 1_i/2., {"Cutoff", 1E-4});
	a = mps.link_c();
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[8]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), s[9]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[7]));
	EXPECT_TRUE(hasindex(mps.Psi_c(), a[9]));
	EXPECT_EQ(mps.Psi_c().inds().size(), 4);
	EXPECT_EQ(mps.get_pos(), 8);
}

TEST_F(MyTest, primesitetest){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	//initial state => ground state of -sx
	myLib::myMPS mps(coeffs);

	myvector<Index> s = mps.site_c();
	myvector<Index> a = mps.link_c();
	mps.primeAll();

	myvector<Index> ps = mps.site_c();
	myvector<Index> pa = mps.link_c();

	for(size_t i=1; i<=s.size(); i++){
		EXPECT_EQ(ps[i], prime(s[i]));
	}
	for(size_t i=1; i<=a.size(); i++){
		EXPECT_EQ(pa[i], prime(a[i]));
	}
}

TEST_F(MyTest, decomposePsitest){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	//initial state => ground state of -sx
	myLib::myMPS mps(coeffs);
	mps.move_pos_to(5);
	mps.decomposePsi();

	myvector<Index> s = mps.site_c();
	myvector<Index> a = mps.link_c();

	EXPECT_TRUE(hasindex(mps.Mref_c(5), a[4]));
	EXPECT_TRUE(hasindex(mps.Mref_c(5), a[5]));
	EXPECT_TRUE(hasindex(mps.Mref_c(5), s[5]));
	EXPECT_EQ(mps.Mref_c(5).inds().size(), 3);

	EXPECT_TRUE(hasindex(mps.Mref_c(6), a[5]));
	EXPECT_TRUE(hasindex(mps.Mref_c(6), a[6]));
	EXPECT_TRUE(hasindex(mps.Mref_c(6), s[6]));
	EXPECT_EQ(mps.Mref_c(6).inds().size(), 3);

	mps.move_pos_to(1);
	mps.decomposePsi();

	a = mps.link_c();

	EXPECT_TRUE(hasindex(mps.Mref_c(1), a[1]));
	EXPECT_TRUE(hasindex(mps.Mref_c(1), s[1]));
	EXPECT_EQ(mps.Mref_c(1).inds().size(), 2);

	EXPECT_TRUE(hasindex(mps.Mref_c(2), a[1]));
	EXPECT_TRUE(hasindex(mps.Mref_c(2), a[2]));
	EXPECT_TRUE(hasindex(mps.Mref_c(2), s[2]));
	EXPECT_EQ(mps.Mref_c(2).inds().size(), 3);

	mps.move_pos_to(9);
	mps.decomposePsi();

	a = mps.link_c();

	EXPECT_TRUE(hasindex(mps.Mref_c(9), a[8]));
	EXPECT_TRUE(hasindex(mps.Mref_c(9), a[9]));
	EXPECT_TRUE(hasindex(mps.Mref_c(9), s[9]));
	EXPECT_EQ(mps.Mref_c(9).inds().size(), 3);

	EXPECT_TRUE(hasindex(mps.Mref_c(10), a[9]));
	EXPECT_TRUE(hasindex(mps.Mref_c(10), s[10]));
	EXPECT_EQ(mps.Mref_c(10).inds().size(), 2);
}

TEST_F(MyTest, timeevolvetest){
	myLib::myvector<std::pair<double, double>> coeffs(10);
	for(auto& coeff : coeffs){
		coeff.first = sqrt(0.5);
		coeff.second = sqrt(0.5);
	}

	//initial state => ground state of -sx
	myLib::myMPS mps(coeffs);

	//initialize each coefficient
	myvector<double> J = std::vector<double>(9, -1);
	myvector<double> h = std::vector<double>(10, -1);
	myvector<double> G = std::vector<double>(10, -1);

	auto A = [](double t){return 1-t;};
	auto B = [](double t){return t;};

	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});
	timeevolve(mps, J, h, G, A, B, 0, 0.1, {"Cutoff",1E-4});

	myvector<ITensor> op_list(10);

	for(size_t i=1; i<=10; i++){
		op_list[i] = id(mps.site_c(i));
	}

	myMPS mps2 = mps;

	EXPECT_NEAR(norm(overlap(mps, op_list, mps2)), 1, 1E-5);

}

TEST_F(MyTest, LZtest){
	double J = 0.75;
	double g = 0.56;
	double T = 100;

	size_t size = 4;

	myLib::myvector<std::pair<double, double>> coeffs_p(size);
	for(auto& coeff : coeffs_p){
		double norm = 1./sqrt(1 +pow((sqrt((g*T)*(g*T) + J*J) + J)/(g*T), 2));
		coeff.first = norm * 1;
		coeff.second = norm * (sqrt((g*T)*(g*T) + J*J) + J)/(g*T);
	}
	myLib::myvector<std::pair<double, double>> coeffs_m(size);
	for(auto& coeff : coeffs_m){
		double norm = 1./sqrt(1 +pow((sqrt((g*T)*(g*T) + J*J) - J)/(g*T), 2));
		coeff.first = norm * 1;
		coeff.second = norm * (sqrt((g*T)*(g*T) + J*J) - J)/(g*T);
	}

	myLib::myMPS mps_p(coeffs_p);
	myLib::myMPS mps_m(coeffs_m, mps_p.site_c());

	//initialize each coefficient
	myvector<double> Jv = std::vector<double>(size-1, 0);
	myvector<double> h = std::vector<double>(size, J/2.);
	myvector<double> G = std::vector<double>(size, g/2.);

	double dt = 0.01;

	auto A = [T](double t){return t;};
	auto B = [T](double t){return 1;};

	for(double t=-T; t<T; t+=dt){
		timeevolve(mps_p, Jv, h, G, A, B, t, dt, {"Cutoff", 1E-4});
	}

	myvector<ITensor> op_list(size);

	for(size_t i=1; i<=size; i++){
		op_list[i] = id(mps_p.site_c(i));
	}

	double Pgs = norm(overlap(mps_m, op_list, mps_p, {"Cutoff", 1E-4})); 
	EXPECT_NEAR(Pgs, exp(-3.14159265*size*(J*J)/(2*g)), 1E-4);
}
