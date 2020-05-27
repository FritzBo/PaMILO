//
//  Point_test.cpp
//  mco
//
//  Created by Fritz Bökler on 30.03.14.
//
//

#include <limits>
#include <tuple>
#include <cmath>

using std::get;

#include <gtest/gtest.h>

using ::testing::Combine;
using ::testing::Values;

#include <mco/basic/point.h>

using mco::Point;

/*********************************************************************
 Construct by Dimension
 ------------------
 Input:
 uint d : dimension
 
 Classes for d:
 (1) 0 < d <= max_int
 (2) d = 0
 
 Expectation
 d          Description
 1          p[i] = 0.0 f.a. i, p.dimension = d
 2          p.values = nullptr, p.dimension = 0
********************************************************************/

TEST(PointTest, ConstructionByDimension) {
    
    unsigned dim = 5;
    Point p1(dim);
    
    EXPECT_EQ(p1.dimension(), dim);
    
    for(unsigned i = 0; i < dim; ++i) {
        EXPECT_EQ(p1[i], 0.0);
    }
    
    Point p2(0);
    
    EXPECT_EQ(p2.dimension(),   0);
    EXPECT_EQ(p2.cbegin(),      nullptr);
}

/*********************************************************************
 Construct by Array
 ------------------
 Two inputs:
 uint d : dimension
 double *a : Data to import
 
 Classes for d:
 (1) 0 < d <= max_int
 (2) d = 0
 
 Classes for a:
 (1) valid chunk of at least d consecutive doubles
 (2) valid chunk of less than d consecutvie doubles
 (3) invalid pointer
 (4) nullptr
 ********************************************************************/

TEST(PointTest, ConstructionByArray) {
    
    // Test for usual input of d: 0 < d < \infty
    unsigned dim = 5;
    
    // -- If a is valid, i.e., points to memory with at least
    // -- d consecutive double values
    double* a = new double[6] {0.3, 1E-6, 43E12, 23, 528, 923.123};
    Point p1(a, dim);
    
    EXPECT_EQ(dim, p1.dimension());
    
    for(unsigned i = 0; i < dim; ++i) {
        EXPECT_EQ(p1[i], a[i]);
    }
    
    // 0 < d < \infty
    // -- If a is nullptr
    // undefined
    
    // 0 < d < \infty
    // -- If len(a) < d
    // undefined
    
    // 0 < d < \infty
    // -- If a is invalid
    // undefined

    // d = 0
    // -- a is valid and len(a) > 0
    Point p3(a, 0);
    
    EXPECT_EQ(p3.dimension(),   0);
    EXPECT_EQ(p3.cbegin(),      nullptr);
    
    // d = 0
    // -- If len(a) < d
    // impossible
    
    // d = 0
    // -- a invalid
    delete[] a;
    Point p4(a, 0);
    
    EXPECT_EQ(0,        p4.dimension());
    EXPECT_EQ(nullptr,  p4.cbegin());

    // d = 0
    // -- a is nullptr
    a = nullptr;
    
    Point p2(a, 0);
    
    EXPECT_EQ(p2.dimension(), 0);
    EXPECT_EQ(p2.cbegin(), nullptr);
}


/*********************************************************************
 Construct by Single Value
 -------------------------
 
 Two inputs:
 uint d : dimension
 double a : value to copy
 
 Classes for d:
 (1) 1 < d <= max_int
 (2) d = 0
 
 Classes for a:
 (1) d is a positive double value
 (2) d is a negative double value
 (3) d is quiet NaN
 (4) d is signaling NaN
 (5) d is infinity
 (6) d is negative infinty
 
 Expectation
 d      a       Description
 1      1       p[i] = a for all i, p.dimension = d
 1      2       same
 1      3       p[i] != a for all i, p.dimension = d
 1      4       same
 1      5       p[i] = a for all i, p.dimension = d
 1      6       same
 
 2      1       p.values = nullptr, p.dimension = 0
 2      2       same
 2      3       same
 2      4       same
 2      5       same
 2      6       same
********************************************************************/

class ConstructBySingleValueFixture
: public ::testing::TestWithParam<std::tuple<unsigned, double>> {
    
public:
    ConstructBySingleValueFixture() : p(test_value, dimension) { }
    
protected:
    unsigned dimension = get<0>(GetParam());
    double test_value = get<1>(GetParam());
    
    Point p;
    
};

TEST_P(ConstructBySingleValueFixture, CheckValue) {
    if(test_value == test_value) {
        for(unsigned i = 0; i < dimension; ++i) {
            EXPECT_EQ(test_value, p[i]);
        }
    } else {
        for(unsigned i = 0; i < dimension; ++i) {
            EXPECT_NE(test_value, p[i]);
        }
    }
}

TEST_P(ConstructBySingleValueFixture, CheckValuesPointer) {
    if(dimension == 0) {
        EXPECT_EQ(nullptr,  p.cbegin());
    } else {
        EXPECT_NE(nullptr, p.cbegin());
    }
}

TEST_P(ConstructBySingleValueFixture, CheckDimension) {
    EXPECT_EQ(dimension, p.dimension());
}

INSTANTIATE_TEST_CASE_P(ConstructBySingleValueTest,
                ConstructBySingleValueFixture,
                Combine(
                        Values(0, 12),
                        Values(3940.23, -1294.232,
                               std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::signaling_NaN(),
                               std::numeric_limits<double>::infinity(),
                               - std::numeric_limits<double>::infinity())
                        )
                );

/*********************************************************************
 Construct by Initializer List
 -----------------------------
 
 One input:
 initializer_list<double> li : initializers
 
 Classes for li:
 (1) li is empty
 (2) lis is a set of double values
 
 Expectation
 li         Description
 1          p.values = nullptr, p.dimension = 0
 2          p[i] = li[i], p.dimension = li.size
 ********************************************************************/
TEST(PointTest, ConstructByInitList) {
    unsigned dim = 5;
    Point p1 {0.3, 1E-6, 43E12, 23, 528};
    
    EXPECT_EQ(dim, p1.dimension());

    EXPECT_EQ(p1[0], 0.3);
    EXPECT_EQ(p1[1], 1E-6);
    EXPECT_EQ(p1[2], 43E12);
    EXPECT_EQ(p1[3], 23);
    EXPECT_EQ(p1[4], 528);
    
    Point p3({});
    
    EXPECT_EQ(p3.dimension(), 0);
    EXPECT_EQ(p3.cbegin(), nullptr);
}

/*********************************************************************
 Copy Construct
 --------------
 
 One input:
 Point p: other point
 
 Classes for p:
 (1) p.dimension = 0
 (2) p.dimension > 0
 
 Expectation
 p          Description
 1          this->values = nullptr, this->dimension = 0
 2          (*this)[i] = p[i], this->dimension = p.dimension
 ********************************************************************/

TEST(PointTest, CopyConstruct) {
    unsigned dim = 5;
    double a = 434.234;
    Point p1(a, 5);
    Point p2(p1);
    
    EXPECT_EQ(dim, p2.dimension());
    EXPECT_NE(p1.cbegin(), p2.cbegin());
    EXPECT_EQ(a, p2[3]);
    
    Point p3(0);
    Point p4(p3);
    
    EXPECT_EQ(0, p4.dimension());
    EXPECT_EQ(nullptr, p4.cbegin());
}

/*********************************************************************
 Move Construct
 --------------
 
 One input:
 Point p: other point
 
 Implicit:
 Point p' after operation
 this pointer
 
 Classes for p:
 (1) p.dimension = 0
 (2) p.dimension > 0
 
 Expectation
 p          Description
 1          this->values = nullptr, this->dimension = 0
 2          (*this)[i] = p[i], this->dimension = p.dimension,
            this->values = p.values, p'.dimension = 0,
            p'. values = nullptr
 ********************************************************************/

TEST(PointTest, MoveConstruct) {
    unsigned dim = 5;
    double a = 434.234;
    
    Point p1(a, 5);
    const double * adress = p1.cbegin();
    
    Point p2(std::move(p1));
    
    EXPECT_EQ(dim,      p2.dimension());
    EXPECT_EQ(adress,   p2.cbegin());
    EXPECT_EQ(a,        p2[3]);
    
    EXPECT_EQ(nullptr,  p1.cbegin());
    EXPECT_EQ(0,        p1.dimension());
    
    Point p3(0);
    Point p4(std::move(p3));
    
    EXPECT_EQ(0, p4.dimension());
    EXPECT_EQ(nullptr, p4.cbegin());
    EXPECT_EQ(nullptr, p4.cbegin());
    EXPECT_EQ(0, p4.dimension());
}

/*********************************************************************
 Copy Assignment
 ---------------
 
 One input:
 Point p: other point
 
 Implicit:
 this pointer
 
 Classes for p:
 (1) p.dimension = 0
 (2) 0 < p.dimension <= max_int
 
 Classes for *this:
 (1) this->dimension = d
 (2) d < this->dimension < max_int
 (3) 0 < this->dimension < d
 (4) this->dimension = 0
 
 Expectation
 p      *this       Description
 1      1           this->values = nullptr, this->dimension = 0
 1      2           this->dimension = 0
 1      3           impossible
 1      4           same as (1-1)
 
 2      1           (*this)[i] = p[i], this->dimension = p.dimension,
                    p.values != this->values
 2      2           same
 2      3           same
 2      4           same
 ********************************************************************/

TEST(PointTest, CopyAssign) {
    // (2-1)
    unsigned dim = 5;
    double a = 434.234;
    Point p1(a, dim);
    Point p2(dim);
    p2 = p1;
    
    EXPECT_EQ(dim, p2.dimension());
    EXPECT_NE(p1.cbegin(), p2.cbegin());
    EXPECT_EQ(a, p2[3]);
    
    // (2-2)
    Point p5(dim - 2);
    p5 = p1;
    EXPECT_EQ(dim, p5.dimension());
    EXPECT_NE(p1.cbegin(), p5.cbegin());
    EXPECT_EQ(a, p5[3]);
    
    // (2-3)
    Point p6(dim + 2);
    p6 = p1;
    EXPECT_EQ(dim, p6.dimension());
    EXPECT_NE(p1.cbegin(), p6.cbegin());
    EXPECT_EQ(a, p6[3]);
    
    // (2-4)
    Point p8(0);
    Point p9(a, dim);
    p8 = p9;
    EXPECT_EQ(dim, p8.dimension());
    EXPECT_NE(p1.cbegin(), p8.cbegin());
    EXPECT_EQ(a, p8[3]);
    
    // (1-1), (1-4)
    Point p3(0);
    Point p4(0);
    p4 = p3;
    
    EXPECT_EQ(0, p4.dimension());
    EXPECT_EQ(nullptr, p4.cbegin());
    
    // (1-2)
    Point p7(3);
    p7 = p3;
    
    EXPECT_EQ(0, p7.dimension());
    
}

/*********************************************************************
 Move Assignment
 ---------------
 
 One input:
 Point p: other point
 
 Implicit:
 this pointer
 
 Classes for p:
 (1) p.dimension = 0
 (2) p.dimension > 0
 
 Expectation
 p          Description
 1          this->values = nullptr, this->dimension = 0
 2          (*this)[i] = p[i], this->dimension = p.dimension,
            p.values != this->values
 ********************************************************************/

TEST(PointTest, MoveAssign) {
    unsigned dim = 5;
    double a = 434.234;
    
    Point p1(a, dim);
    const double * adress = p1.cbegin();
    
    Point p2(dim);
    p2 = std::move(p1);
    
    EXPECT_EQ(dim, p2.dimension());
    EXPECT_EQ(adress, p2.cbegin());
    EXPECT_EQ(a, p2[3]);
    
    Point p3(0);
    Point p4(0);
    p4 = std::move(p3);
    
    EXPECT_EQ(0, p4.dimension());
    EXPECT_EQ(nullptr, p4.cbegin());
    EXPECT_EQ(nullptr, p4.cbegin());
    EXPECT_EQ(0, p4.dimension());
}

TEST(PointTest, AddInplace) {
    Point p1 = {1, 2, 3};
    Point p2 = {3, 2, 1};
    
    p1 += p2;
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(4, p1[i]);
        EXPECT_EQ(3 - i, p2[i]);
    }
    
    Point p3(0);
    Point p4(0);
    
    p3 += p4;
    
    EXPECT_EQ(0, p3.dimension());
    EXPECT_EQ(nullptr, p3.cbegin());
    
#ifdef DEBUG
    
    EXPECT_DEBUG_DEATH(p1 += p3, "Assertion failed");
    EXPECT_DEBUG_DEATH(p3 += p1, "Assertion failed");
    
#endif
}

TEST(PointTest, SubtractInplace) {
    Point p1 = {1, 2, 3};
    Point p2 = {3, 2, 1};
    
    p1 -= p2;
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(-2 + 2 * i,   p1[i]);
        EXPECT_EQ(3 - i,        p2[i]);
    }
    
    Point p3(0);
    Point p4(0);
    
    p3 -= p4;
    
    EXPECT_EQ(0,        p3.dimension());
    EXPECT_EQ(nullptr,  p3.cbegin());
    
#ifdef DEBUG
    
    EXPECT_DEBUG_DEATH(p1 -= p3, "Assertion failed");
    EXPECT_DEBUG_DEATH(p3 -= p1, "Assertion failed");
    
#endif

    
}

TEST(PointTest, ScalarProdInplace) {
    double d1 = 5;
    double d2 = -10;
    
    Point p = {1,2,3,4,5,6};
    
    p *= d1;
    
    for(auto i : {0, 1, 2, 3, 4, 5}) {
        EXPECT_EQ((i + 1) * 5,  p[i]);
    }
    
    p *= d2;
    
    for(auto i : {0, 1, 2, 3, 4, 5}) {
        EXPECT_EQ((i + 1) * 5 * -10,  p[i]);
    }
    
    Point p2(0);
    
    p2 *= d1;
    
    EXPECT_EQ(0,        p2.dimension());
    EXPECT_EQ(nullptr,  p2.cbegin());
}

TEST(PointTest, Negation) {
    Point p1 = {1, 2, 3, 4, 5};
    
    for(auto i : {0, 1, 2, 3, 4}) {
        EXPECT_EQ(-p1[i], (-p1)[i]);
    }
    
    for(auto i : {0, 1, 2, 3, 4}) {
        EXPECT_EQ(i + 1, p1[i]);
    }
    
    const Point p2(0);
    
    EXPECT_EQ(0,        (-p2).dimension());
    EXPECT_EQ(nullptr,  (-p2).cbegin());
}

TEST(PointTest, NegationRValue) {
    Point p = {1, 2, 3, 4, 5};
    
    Point p1(p);
    const double * adress = p1.cbegin();
    
    Point p2 = -std::move(p1);
    
    EXPECT_EQ(adress, p2.cbegin());
    
    for(auto i : {0, 1, 2, 3, 4}) {
        EXPECT_EQ(-p[i], p2[i]);
    }
    
    Point p3(0);
    p3 = - Point(0);
    
    EXPECT_EQ(0,        p3.dimension());
    EXPECT_EQ(nullptr,  p3.cbegin());
}

TEST(PointTest, Addition) {
    unsigned dim = 3;
    
    Point p1 = {134.2, 23, 223.4};
    Point p2 = {1.9, 343.2, 123.34};
    
    Point p3(p2);
    const double * p3_address = p3.cbegin();
    
    Point p4(p1);
    const double * p4_address = p4.cbegin();
    
    Point p5(p1);
    const double * p5_address = p5.cbegin();
    
    Point p6(p2);
    const double * p6_address = p6.cbegin();
    
    Point p7    = p1            + p2;
    Point p8    = p1            + std::move(p3);
    Point p9    = std::move(p4) + p2;
    Point p10   = std::move(p5) + std::move(p6);
    
    EXPECT_EQ(dim, p7.dimension());
    EXPECT_EQ(dim, p8.dimension());
    EXPECT_EQ(dim, p9.dimension());
    EXPECT_EQ(dim, p10.dimension());
    
    EXPECT_NE(p1.cbegin(), p7.cbegin());
    EXPECT_NE(p2.cbegin(), p7.cbegin());
    
    EXPECT_NE(p1.cbegin(), p8.cbegin());
    EXPECT_EQ(p3_address,   p8.cbegin());
    
    EXPECT_NE(p2.cbegin(), p9.cbegin());
    EXPECT_EQ(p4_address,   p9.cbegin());
    
    EXPECT_TRUE(p10.cbegin() == p5_address ||
                p10.cbegin() == p6_address);
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(p1[i] + p2[i], p7[i]);
    }
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(p1[i] + p2[i], p8[i]);
    }
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(p1[i] + p2[i], p9[i]);
    }
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(p1[i] + p2[i], p10[i]);
    }
    
#ifdef DEBUG
    
    Point p11(0);
    
    EXPECT_DEBUG_DEATH(p1       + p11,      "Assertion failed");
    EXPECT_DEBUG_DEATH(p1       + Point(6), "Assertion failed");
    EXPECT_DEBUG_DEATH(Point(2) + p1,       "Assertion failed");
    EXPECT_DEBUG_DEATH(Point(2) + Point(8), "Assertion failed");
    
#endif
}

TEST(PointTest, Subtraction) {
    unsigned dim = 3;
    
    Point p1 = {134.2, 23, 223.4};
    Point p2 = {1.9, 343.2, 123.34};
    
    Point p3(p2);
    const double * p3_address = p3.cbegin();
    
    Point p4(p1);
    const double * p4_address = p4.cbegin();
    
    Point p5(p1);
    const double * p5_address = p5.cbegin();
    
    Point p6(p2);
    const double * p6_address = p6.cbegin();
    
    Point p7 = p1 - p2;
    Point p8 = p1 - std::move(p3);
    Point p9 = std::move(p4) - p2;
    Point p10 = std::move(p5) - std::move(p6);
    
    EXPECT_EQ(dim, p7.dimension());
    EXPECT_EQ(dim, p8.dimension());
    EXPECT_EQ(dim, p9.dimension());
    EXPECT_EQ(dim, p10.dimension());
    
    EXPECT_NE(p1.cbegin(), p7.cbegin());
    EXPECT_NE(p2.cbegin(), p7.cbegin());
    
    EXPECT_NE(p1.cbegin(), p8.cbegin());
    EXPECT_EQ(p3_address, p8.cbegin());
    
    EXPECT_NE(p2.cbegin(), p9.cbegin());
    EXPECT_EQ(p4_address, p9.cbegin());
    
    EXPECT_TRUE(p10.cbegin() == p5_address ||
                p10.cbegin() == p6_address);
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(p1[i] - p2[i], p7[i]);
    }
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(p1[i] - p2[i], p8[i]);
    }
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(p1[i] - p2[i], p9[i]);
    }
    
    for(auto i : {0, 1, 2}) {
        EXPECT_EQ(p1[i] - p2[i], p10[i]);
    }
    
#ifdef DEBUG
    
    Point p11(0);
    
    EXPECT_DEBUG_DEATH(p1       - p11,      "Assertion failed");
    EXPECT_DEBUG_DEATH(p1       - Point(6), "Assertion failed");
    EXPECT_DEBUG_DEATH(Point(2) - p1,       "Assertion failed");
    EXPECT_DEBUG_DEATH(Point(2) - Point(8), "Assertion failed");
    
#endif

}

TEST(PointTest, ScalarProduct) {
    double d1 = 543.234;
    double d2 = -10;
    
    Point p1 = {594.309, 92.32, 1.32523, 3.42};
    
    Point p2 = p1 * d1;
    
    for(unsigned i = 0; i < p1.dimension(); ++i) {
        EXPECT_EQ(p1[i] * d1,  p2[i]);
    }
    
    Point p3 = p1 * d2;
    
    for(unsigned i = 0; i < p1.dimension(); ++i) {
        EXPECT_EQ(p1[i] * d2,  p3[i]);
    }
    
    Point p4(0);
    
    Point p5 = p4 * d1;
    
    EXPECT_EQ(0,        p5.dimension());
    EXPECT_EQ(nullptr,  p5.cbegin());

}

TEST(PointTest, InnerProduct) {
    Point p1 = {0, 0, 1};
    Point p2 = {1, 0, 0};
    
    EXPECT_EQ(0, p1 * p2);
    
    p1 = {123, 1.32, 229.23};
    p2 = {12.4, 123.93, 934.1};
    
    EXPECT_EQ(123 * 12.4 +
              1.32 * 123.93 +
              229.23 * 934.1, p1 * p2);
    
#ifdef DEBUG
    EXPECT_DEBUG_DEATH(p1 * Point(5), "Assertion failed");
#endif
    
}

TEST(PointTest, OperatorLeftShift) {
    EXPECT_NO_FATAL_FAILURE(std::cout << Point(0) << std::endl);
    EXPECT_NO_FATAL_FAILURE(std::cout << Point({1,2,3}) << std::endl);
    EXPECT_NO_FATAL_FAILURE(std::cout << Point({8953.12, 835.3, 284.84}) << std::endl);
}

TEST(PointTest, ComponentwisePointComparator) {
    double epsilon = 1E-5;
    
    // Base
    Point p1 = {9.2952,     8.28232,    9.92242,    2.1482,     3.22995};
    
    // larger
    Point p2 = {10.2529,    9.12529,    12.29523,   3.1253,     3.91235};
    
    // larger and equal
    Point p3 = {12.1295,    8.28232,    11.49252,   4.2952,     5.29523};
    
    // larger and smaller
    Point p4 = {11.1295,    5.29523,    13.12492,   2.9124,     4.52343};
    
    // larger and not large enough
    Point p5 = {p1[0] + epsilon / 2,
                            9.23942,    12.1239,    5.1239,     8.12394};
    
    // larger and large enough
    Point p6 = {15.23923,   9.3123,     p1[2] + 2 * epsilon,
                                                    4.2382,     4.23923};
    
    mco::ComponentwisePointComparator comp_le;
    
    EXPECT_FALSE(comp_le(p1, p1));
    EXPECT_TRUE(comp_le(p1, p2));
    EXPECT_FALSE(comp_le(p1, p3));
    EXPECT_FALSE(comp_le(p1, p4));
    EXPECT_TRUE(comp_le(p1, p5));
    EXPECT_TRUE(comp_le(p1, p6));
    
    mco::ComponentwisePointComparator comp_e_le(epsilon);
    
    EXPECT_FALSE(comp_e_le(p1, p1));
    EXPECT_TRUE(comp_e_le(p1, p2));
    EXPECT_FALSE(comp_e_le(p1, p3));
    EXPECT_FALSE(comp_e_le(p1, p4));
    EXPECT_FALSE(comp_e_le(p1, p5));
    EXPECT_TRUE(comp_e_le(p1, p6));
    
    mco::ComponentwisePointComparator comp_leq(0, false);
    
    EXPECT_TRUE(comp_leq(p1, p1));
    EXPECT_TRUE(comp_leq(p1, p2));
    EXPECT_TRUE(comp_leq(p1, p3));
    EXPECT_FALSE(comp_leq(p1, p4));
    EXPECT_TRUE(comp_leq(p1, p5));
    EXPECT_TRUE(comp_leq(p1, p6));
    
    mco::ComponentwisePointComparator comp_e_leq(epsilon, false);
    
    EXPECT_TRUE(comp_e_leq(p1, p1));
    EXPECT_TRUE(comp_e_leq(p1, p2));
    EXPECT_TRUE(comp_e_leq(p1, p3));
    EXPECT_FALSE(comp_e_leq(p1, p4));
    EXPECT_TRUE(comp_e_leq(p1, p5));
    EXPECT_TRUE(comp_e_leq(p1, p6));
}

TEST(PointTest, EqualityPointComparator) {
    double epsilon = 1E-5;
    
    // base
    Point p1 = {9.2952,     8.28232,    9.92242,    2.1482,     3.22995};
    
    // different
    Point p2 = {2.5233,     9.32832,    9.32134,    4.1293,     4.23942};
    
    // some equal
    Point p3 = {8.3293,     p1[1],      5.23942,    p1[3],      9.23923};
    
    // almost same
    Point p5 = {p1[0] + epsilon / 2,
                p1[1] - epsilon / 3,
                p1[2] + epsilon / 5,
                p1[3] - epsilon / 10,
                p1[4] + epsilon / 7};
    
    // almost same but different
    Point p6 = {p1[0] + epsilon / 3,
                p1[1] - epsilon / 2,
                p1[2] + epsilon / 5,
                p1[3] - epsilon * 2,
                p1[4] + epsilon / 7};
    
    mco::EqualityPointComparator eq;
    
    EXPECT_TRUE(eq(p1, p1));
    EXPECT_FALSE(eq(p1, p2));
    EXPECT_FALSE(eq(p1, p3));
    EXPECT_TRUE(eq(p1, p1));
    EXPECT_FALSE(eq(p1, p5));
    EXPECT_FALSE(eq(p1, p6));
    
    mco::EqualityPointComparator eq_e(epsilon);
    
    EXPECT_TRUE(eq_e(p1, p1));
    EXPECT_FALSE(eq_e(p1, p2));
    EXPECT_FALSE(eq_e(p1, p3));
    EXPECT_TRUE(eq_e(p1, p1));
    EXPECT_TRUE(eq_e(p1, p5));
    EXPECT_FALSE(eq_e(p1, p6));
    
}

TEST(PointTest, LexPointComparator) {
    double epsilon = 1E-5;
    
    // base
    Point p1 = {9.2952,     8.28232,    9.92242,    2.1482,     3.22995};
    
    // larger in all components
    Point p2 = {11.293,     9.23923,    10.2392,    5.2392,     4.23923};
    
    // larger in first component, smaller in all other
    Point p3 = {12.2392,    2.29323,    5.23932,    1.2392,     2.3233};
    
    // smaller in first compoenent, larger in all other
    Point p4 = {5.12392,    12.1239,    10.2392,    5.12392,    8.1239};
    
    // same in first two components, larger in third, smaller in rest
    Point p5 = {p1[0],      p1[1],      10.2932,    1.9232,     1.9232};
    
    // slightly smaller in first two components, larger in third
    Point p6 = {p1[0] - epsilon / 2,
                p1[1] - epsilon / 5,
                p1[2] + epsilon * 10,
                p1[3] - epsilon * 42,
                p1[4] - epsilon * 534};
    
    // almost same but lex-larger
    Point p7 = {p1[0] + epsilon / 2,
                p1[1] - epsilon / 3,
                p1[2] + epsilon / 5,
                p1[3] - epsilon / 10,
                p1[4] + epsilon / 7};
    
    // almost same but lex-smaller
    Point p8 = {p1[0] - epsilon / 2,
                p1[1] + epsilon / 3,
                p1[2] - epsilon / 5,
                p1[3] + epsilon / 10,
                p1[4] - epsilon / 7};
    
    // almost same in first two components, larger in third, smaller in rest
    Point p9 = {p1[0] - epsilon / 2,
                p1[1] - epsilon / 3,
                                        10.2932,    1.9232,     1.9232};

    
    mco::LexPointComparator lex_le;
    
    EXPECT_FALSE(lex_le(p1, p1));
    EXPECT_TRUE(lex_le(p1, p2));
    EXPECT_TRUE(lex_le(p1, p3));
    EXPECT_FALSE(lex_le(p1, p4));
    EXPECT_TRUE(lex_le(p1, p5));
    EXPECT_FALSE(lex_le(p1, p6));
    EXPECT_TRUE(lex_le(p1, p7));
    EXPECT_FALSE(lex_le(p1, p8));
    EXPECT_FALSE(lex_le(p1, p9));
    
    mco::LexPointComparator lex_e_le(epsilon);
    
    EXPECT_FALSE(lex_e_le(p1, p1));
    EXPECT_TRUE(lex_e_le(p1, p2));
    EXPECT_TRUE(lex_e_le(p1, p3));
    EXPECT_FALSE(lex_e_le(p1, p4));
    EXPECT_TRUE(lex_e_le(p1, p5));
    EXPECT_TRUE(lex_e_le(p1, p6));
    EXPECT_FALSE(lex_e_le(p1, p7));
    EXPECT_FALSE(lex_e_le(p1, p8));
    EXPECT_TRUE(lex_e_le(p1, p9));
    
    mco::LexPointComparator lex_leq(0, false);
    
    EXPECT_TRUE(lex_leq(p1, p1));
    EXPECT_TRUE(lex_leq(p1, p2));
    EXPECT_TRUE(lex_leq(p1, p3));
    EXPECT_FALSE(lex_leq(p1, p4));
    EXPECT_TRUE(lex_leq(p1, p5));
    EXPECT_FALSE(lex_leq(p1, p6));
    EXPECT_TRUE(lex_le(p1, p7));
    EXPECT_FALSE(lex_le(p1, p8));
    EXPECT_FALSE(lex_le(p1, p9));
    
    mco::LexPointComparator lex_e_leq(epsilon, false);
    
    EXPECT_TRUE(lex_e_leq(p1, p1));
    EXPECT_TRUE(lex_e_leq(p1, p2));
    EXPECT_TRUE(lex_e_leq(p1, p3));
    EXPECT_FALSE(lex_e_leq(p1, p4));
    EXPECT_TRUE(lex_e_leq(p1, p5));
    EXPECT_TRUE(lex_e_leq(p1, p6));
    EXPECT_TRUE(lex_e_leq(p1, p7));
    EXPECT_TRUE(lex_e_leq(p1, p8));
    EXPECT_TRUE(lex_e_leq(p1, p9));
}


class ParetoPointComparatorFixture : public ::testing::TestWithParam<std::tuple<Point, Point>> {
    
};

TEST_P(ParetoPointComparatorFixture, Dominates) {
    
}

TEST_P(ParetoPointComparatorFixture, NonDominated) {
    
}

TEST_P(ParetoPointComparatorFixture, Incomparable) {
    
}

INSTANTIATE_TEST_CASE_P(ParetoPointComparatorTest, ParetoPointComparatorFixture, Combine(
                        Values(Point(0), Point(5), Point(10)),
                        Values(Point(0), Point(3), Point(5))));

TEST(PointTest, ParetoPointComparator) {
    double epsilon = 1E-5;
    
    // base
    Point p = {9.2952,     8.28232,    9.92242,    2.1482,     3.22995};
    
    // dominated, non-weak
    Point p1 = {p[0] + epsilon * 2323,
                p[1] + epsilon * 423,
                p[2] + epsilon * 9423,
                p[3] + epsilon * 5823,
                p[4] + epsilon * 4293};
    
    // dominated, weak, exact
    Point p2 = {p[0] + epsilon * 2323,
                p[1] + epsilon * 423,
                p[2],
                p[3] + epsilon * 5823,
                p[4] + epsilon * 4293};

    
    // dominated, almost same
    Point p3 = {p[0] + epsilon / 2,
                p[1] + epsilon / 4,
                p[2] + epsilon / 2,
                p[3] + epsilon / 2,
                p[4] + epsilon / 4};
    
    // dominated, some same
    Point p4 = {p[0] + epsilon / 2,
                p[1] + epsilon / 4,
                p[2] + epsilon * 5329,
                p[3] + epsilon / 2,
                p[4] + epsilon * 5394};


    
    // incomparable, no epsilon
/*    Point p4 = {p1[0] + epsilon * 9233,
                p1[1] - epsilon * 2932,
                p1[2] - epsilon * 4293,
                p1[3] - epsilon * 9423,
                p1[4] + epsilon * 4232};
    
    // non-dominated, no epsilon
    Point p5 = {p1[0] - epsilon * 2323,
                p1[1] - epsilon * 423,
                p1[2] - epsilon * 9423,
                p1[3] - epsilon * 5823,
                p1[4] - epsilon * 4293};
    
    // incomparable, almost same
    Point p8 = {p1[0] + epsilon / 2,
                p1[1] + epsilon / 4,
                p1[2] + epsilon / 10,
                p1[3] - epsilon / 8,
                p1[4] - epsilon / 3};
    
    // incomparable, some same
    Point p9 = {p1[0] - epsilon / 2,
                p1[1] + epsilon / 4,
                p1[2] - epsilon / 2,
                p1[3] + epsilon / 2,
                p1[4] + epsilon / 4};*/

    
    mco::ParetoDominationPointComparator pareto_le;
    
    
    
    mco::ParetoDominationPointComparator pareto_le_e(epsilon);
}