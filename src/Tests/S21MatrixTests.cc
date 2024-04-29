#include <fcntl.h>
#include <gtest/gtest.h>
#include <cmath>
#include "../Lib/S21Matrix.h"


void fillMatrix(S21Matrix& matrix) {
  for (int i = 0; i < matrix.getRows(); ++i) {
    for (int j = 0; j < matrix.getCols(); ++j) {
      matrix[i][j] = rand() % 10;
    }
  }
}

void setIdentityMatrix(S21Matrix& matrix) {
  for (int i = 0; i < matrix.getRows(); ++i) {
    for (int j = 0; j < matrix.getCols(); ++j) {
      matrix[i][j] = !(i ^ j);
    }
  }
}

TEST(EDefaultConstructorExceptionTest, test1) {
  EXPECT_ANY_THROW(S21Matrix(3, 0));
  EXPECT_ANY_THROW(S21Matrix(3, -1));
  EXPECT_ANY_THROW(S21Matrix(-1, 3));
  EXPECT_ANY_THROW(S21Matrix(0, 3));
  EXPECT_ANY_THROW(S21Matrix(0, 0));
  EXPECT_ANY_THROW(S21Matrix(-1, -1));
}

TEST(OperatorsExceptionTest, test1) {
  S21Matrix a(1, 1);
  S21Matrix b(2, 2);

  EXPECT_ANY_THROW(a + b);
  EXPECT_ANY_THROW(a - b);
  EXPECT_ANY_THROW(a * b);
  EXPECT_ANY_THROW(a += b);
  EXPECT_ANY_THROW(a -= b);
  EXPECT_ANY_THROW(a *= b);
}

TEST(MatrixOperationsExceptionTest, test1) {
  S21Matrix a(1, 2);
  S21Matrix c(3, 3);

  for (int i = 0; i < c.getRows(); ++i) {
    for (int j = 0; j < c.getCols(); ++j) {
      c[i][j] = 1;
    }
  }

  EXPECT_ANY_THROW(a.CalcComplements());
  EXPECT_ANY_THROW(a.InverseMatrix());
  EXPECT_ANY_THROW(c.CalcComplements());
  EXPECT_ANY_THROW(c.InverseMatrix());
}

TEST(GetterTest, test1) {
  S21Matrix a(3, 3);
  EXPECT_EQ(a.getRows(), 3);
  EXPECT_EQ(a.getCols(), 3);

  S21Matrix b(1, 3);
  EXPECT_EQ(b.getRows(), 1);
  EXPECT_EQ(b.getCols(), 3);

  b.resize(1, 1);
  EXPECT_EQ(b.getRows(), 1);
  EXPECT_EQ(b.getCols(), 1);
}

TEST(SetterTest, test1) {
  S21Matrix a(3, 3);
  fillMatrix(a);
  S21Matrix b(a);

  a.resize(10, 10);
  a.resize(3, 3);
  a.resize(3, 3);

  EXPECT_TRUE(a == b);
}

TEST(SetterTest, test2) {
  S21Matrix a(3, 3);
  a[0][0] = 1;
  a[0][1] = 2;
  a[0][2] = 3;

  a.setCols(5);
  EXPECT_EQ(a.getCols(), 5);

  a.setRows(5);
  EXPECT_EQ(a.getRows(), 5);

  EXPECT_DOUBLE_EQ(a[0][0], 1);
  EXPECT_DOUBLE_EQ(a[0][1], 2);
  EXPECT_DOUBLE_EQ(a[0][2], 3);

  a.resize(1, 1);
  EXPECT_EQ(a.getCols(), 1);
  EXPECT_EQ(a.getRows(), 1);
}

TEST(BracketOperatorTest, test1) {
  S21Matrix a(3, 3);
  a[0][0] = 1;
  a[0][1] = 2;
  a[0][2] = 3;
  a[1][0] = 4;
  a[1][1] = 5;
  a[1][2] = 6;
  a[2][0] = 7;
  a[2][1] = 8;  
  a[2][2] = 9;

  double val = 0;
  for (int i = 0; i < a.getRows(); ++i) {
    for (int j = 0; j < a.getCols(); ++j) {
      EXPECT_DOUBLE_EQ(a[i][j], ++val);
    }
  }
}

TEST(EqualOperatorTest, test1) {
  S21Matrix a(3, 3);
  fillMatrix(a);
  S21Matrix b(a);

  EXPECT_TRUE(a == b);
  b[0][0] = !b[0][0];
  EXPECT_TRUE(a != b);
}

TEST(EqualOperatorTest, test2) {
  S21Matrix a(3, 3);
  fillMatrix(a);
  S21Matrix b(3, 3);

  EXPECT_TRUE(a != b);
}

TEST(ConstructorTest, test1) {
  S21Matrix a;
  ASSERT_EQ(a.getRows(), 0);
  ASSERT_EQ(a.getCols(), 0);
}

TEST(ConstructorTest, test2) {
  S21Matrix a(3, 3);
  ASSERT_EQ(a.getRows(), 3);
  ASSERT_EQ(a.getCols(), 3);
}

TEST(ConstructorTest, test3) {
  S21Matrix a(3, 3);
  fillMatrix(a);
  S21Matrix b(a);
  EXPECT_TRUE(a == b);
}

TEST(ConstructorTest, test4) {
  S21Matrix a(3, 3);
  fillMatrix(a);

  S21Matrix b(a);
  S21Matrix c(std::move(a));

  EXPECT_TRUE(c == b);

  S21Matrix empty;
  EXPECT_TRUE(a == empty);
}

TEST(AssignmentOperatorTest, test1) {
  S21Matrix a(3, 3);
  fillMatrix(a);

  S21Matrix b(3, 5);
  b = a;
  EXPECT_TRUE(a == b);
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(AdditionAssignmentOperatorTest, test1) {
  S21Matrix a(3, 3);
  fillMatrix(a);

  S21Matrix b(3, 3);
  fillMatrix(a);

  S21Matrix result = a;
  result += b;

  EXPECT_TRUE(result == a + b);
}

TEST(SubtractionAssignmentOperatorTest, test1) {
  S21Matrix a(3, 3);
  fillMatrix(a);

  S21Matrix b(3, 3);
  fillMatrix(a);

  S21Matrix result = a;
  result -= b;

  EXPECT_TRUE(result == a - b);
}

TEST(MultiplicationAssignmentOperatorTest, test1) {
  S21Matrix a(3, 3);
  fillMatrix(a);

  S21Matrix b(3, 3);
  fillMatrix(a);

  S21Matrix result = a;
  result *= b;

  EXPECT_TRUE(result == a * b);
}

TEST(MultiplicationAssignmentOperatorTest, test2) {
  S21Matrix a(3, 3);
  fillMatrix(a);

  double value = 1.5;

  S21Matrix result = a;
  result *= value;

  EXPECT_TRUE(result == a * value);
}

TEST(SumTest, test0) {
  S21Matrix a1(3, 3);
  fillMatrix(a1);

  S21Matrix b1(3, 3);
  fillMatrix(b1);

  S21Matrix a2(a1);
  S21Matrix b2(b1);

  a1 += b1;
  a2.SumMatrix(b2);
  EXPECT_TRUE(a1 == a2);
}

TEST(SubTest, test0) {
  S21Matrix a1(3, 3);
  fillMatrix(a1);

  S21Matrix b1(3, 3);
  fillMatrix(b1);

  S21Matrix a2(a1);
  S21Matrix b2(b1);

  a1 -= b1;
  a2.SubMatrix(b2);
  EXPECT_TRUE(a1 == a2);
}

TEST(MulTest, test0) {
  S21Matrix a1(3, 3);
  fillMatrix(a1);

  S21Matrix b1(3, 3);
  fillMatrix(b1);

  S21Matrix a2(a1);
  S21Matrix b2(b1);

  a1 *= b1;
  a2.MulMatrix(b2);
  EXPECT_TRUE(a1 == a2);
}

TEST(MulTest, test2) {
  S21Matrix a1(3, 3);
  fillMatrix(a1);

  double value = 1.23;

  S21Matrix a2(a1);

  a1 *= value;
  a2.MulNumber(value);
  EXPECT_TRUE(a1 == a2);
}

TEST(TransposeTest, test1) {
  S21Matrix a(3, 5);
  fillMatrix(a);

  S21Matrix b(a);

  b = b.Transpose().Transpose();
  EXPECT_TRUE(a == b);
}

TEST(InverseTest, test1) {
    S21Matrix a(3, 3);
    fillMatrix(a);

    EXPECT_TRUE(a == a.InverseMatrix().InverseMatrix());
}

TEST(MinorTest, test1) {
  S21Matrix a(3, 3);
  a[0][0] = 1;
  a[0][1] = 2;
  a[0][2] = 3;
  a[1][0] = 4;
  a[1][1] = 5;
  a[1][2] = 6;
  a[2][0] = 7;
  a[2][1] = 8;
  a[2][2] = 9;

  S21Matrix b = a.Minor(1, 1);

  EXPECT_DOUBLE_EQ(b[0][0], a[0][0]);
  EXPECT_DOUBLE_EQ(b[0][1], a[0][2]);
  EXPECT_DOUBLE_EQ(b[1][0], a[2][0]);
  EXPECT_DOUBLE_EQ(b[1][1], a[2][2]);
}

TEST(DeterminantTest, test1) {
  S21Matrix a(3, 3);
  a[0][0] = 1;
  a[0][1] = 2;
  a[0][2] = 3;
  a[1][0] = 4;
  a[1][1] = 5;
  a[1][2] = 6;
  a[2][0] = 7;
  a[2][1] = 8;
  a[2][2] = 9;

  EXPECT_EQ(a.Determinant(), 0);
}

TEST(DeterminantTest, test2) {
  S21Matrix a(3, 3);
  setIdentityMatrix(a);

  EXPECT_DOUBLE_EQ(a.Determinant(), 1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}