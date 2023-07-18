#include <Eigen/Eigen>
#include <cmath>
#include <vector>

// 与实际的 exp() 和 log() 函数相比，这里的 expC2() 和 logC2() 函数
// 曲线精度差别巨大 不明白为啥要这样用 时间感觉差了一倍
static double expC2(double t) {
  return t > 0.0 ? ((0.5 * t + 1.0) * t + 1.0)
                 : 1.0 / ((0.5 * t - 1.0) * t + 1.0);
}
static double logC2(double T) {
  return T > 1.0 ? (sqrt(2.0 * T - 1.0) - 1.0) : (1.0 - sqrt(2.0 / T - 1.0));
}


// The banded system class is used for solving
// banded linear system Ax=b efficiently.
// A is an N*N band matrix with lower band width lowerBw
// and upper band width upperBw.
// Banded LU factorization has O(N) time complexity.
class BandedSystem {
 public:
  // The size of A, as well as the lower/upper
  // banded width p/q are needed
  inline void create(const int &n, const int &p, const int &q) {
    // In case of re-creating before destroying
    destroy();
    N = n;
    lowerBw = p;
    upperBw = q;
    int actualSize = N * (lowerBw + upperBw + 1);
    ptrData = new double[actualSize];
    std::fill_n(ptrData, actualSize, 0.0);
    return;
  }

  inline void destroy() {
    if (ptrData != nullptr) {
      delete[] ptrData;
      ptrData = nullptr;
    }
    return;
  }

 private:
  int N;
  int lowerBw;
  int upperBw;
  // Compulsory nullptr initialization here
  double *ptrData = nullptr;

 public:
  // Reset the matrix to zero
  inline void reset(void) {
    std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
    return;
  }

  // The band matrix is stored as suggested in "Matrix Computation"
  inline const double &operator()(const int &i, const int &j) const {
    return ptrData[(i - j + upperBw) * N + j];
  }

  inline double &operator()(const int &i, const int &j) {
    return ptrData[(i - j + upperBw) * N + j];
  }

  // This function conducts banded LU factorization in place
  // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
  inline void factorizeLU() {
    int iM, jM;
    double cVl;
    for (int k = 0; k <= N - 2; k++) {
      iM = std::min(k + lowerBw, N - 1);
      cVl = operator()(k, k);
      for (int i = k + 1; i <= iM; i++) {
        if (operator()(i, k) != 0.0) {
          operator()(i, k) /= cVl;
        }
      }
      jM = std::min(k + upperBw, N - 1);
      for (int j = k + 1; j <= jM; j++) {
        cVl = operator()(k, j);
        if (cVl != 0.0) {
          for (int i = k + 1; i <= iM; i++) {
            if (operator()(i, k) != 0.0) {
              operator()(i, j) -= operator()(i, k) * cVl;
            }
          }
        }
      }
    }
    return;
  }

  // This function solves Ax=b, then stores x in b
  // The input b is required to be N*m, i.e.,
  // m vectors to be solved.
  inline void solve(Eigen::MatrixXd &b) const {
    int iM;
    for (int j = 0; j <= N - 1; j++) {
      iM = std::min(j + lowerBw, N - 1);
      for (int i = j + 1; i <= iM; i++) {
        if (operator()(i, j) != 0.0) {
          b.row(i) -= operator()(i, j) * b.row(j);
        }
      }
    }
    for (int j = N - 1; j >= 0; j--) {
      b.row(j) /= operator()(j, j);
      iM = std::max(0, j - upperBw);
      for (int i = iM; i <= j - 1; i++) {
        if (operator()(i, j) != 0.0) {
          b.row(i) -= operator()(i, j) * b.row(j);
        }
      }
    }
    return;
  }

  // This function solves ATx=b, then stores x in b
  // The input b is required to be N*m, i.e.,
  // m vectors to be solved.
  inline void solveAdj(Eigen::MatrixXd &b) const {
    int iM;
    for (int j = 0; j <= N - 1; j++) {
      b.row(j) /= operator()(j, j);
      iM = std::min(j + upperBw, N - 1);
      for (int i = j + 1; i <= iM; i++) {
        if (operator()(j, i) != 0.0) {
          b.row(i) -= operator()(j, i) * b.row(j);
        }
      }
    }
    for (int j = N - 1; j >= 0; j--) {
      iM = std::max(0, j - lowerBw);
      for (int i = iM; i <= j - 1; i++) {
        if (operator()(j, i) != 0.0) {
          b.row(i) -= operator()(j, i) * b.row(j);
        }
      }
    }
    return;
  }
};



// MATLAB code for the function
classdef BandedSystem
    properties
        N
        lowerBw
        upperBw
        ptrData
    end
    
    methods
        function obj = BandedSystem(n, p, q)
            obj.N = n;
            obj.lowerBw = p;
            obj.upperBw = q;
            actualSize = obj.N * (obj.lowerBw + obj.upperBw + 1);
            obj.ptrData = zeros(1, actualSize);
        end
        
        function delete(obj)
            % No need to implement explicit destructor in MATLAB
        end
        
        function reset(obj)
            obj.ptrData = zeros(1, obj.N * (obj.lowerBw + obj.upperBw + 1));
        end
        
        function value = subsref(obj, s)
            % Overload MATLAB's subsref to enable accessing elements using parenthesis
            if strcmp(s(1).type, '()')
                i = s(1).subs{1};
                j = s(1).subs{2};
                value = obj.ptrData((i - j + obj.upperBw) * obj.N + j);
            end
        end
        
        function obj = subsasgn(obj, s, value)
            % Overload MATLAB's subsasgn to enable assigning values using parenthesis
            if strcmp(s(1).type, '()')
                i = s(1).subs{1};
                j = s(1).subs{2};
                obj.ptrData((i - j + obj.upperBw) * obj.N + j) = value;
            end
        end

        function obj = factorizeLU(obj)
            iM = min(k + obj.lowerBw, obj.N - 1);
            for k = 0:(obj.N - 2)
                cVl = obj(k + 1, k + 1);
                for i = (k + 2):iM
                    if obj(i, k + 1) ~= 0.0
                        obj(i, k + 1) = obj(i, k + 1) / cVl;
                    end
                end
                jM = min(k + obj.upperBw, obj.N - 1);
                for j = (k + 2):jM
                    cVl = obj(k + 1, j);
                    if cVl ~= 0.0
                        for i = (k + 2):iM
                            if obj(i, k + 1) ~= 0.0
                                obj(i, j) = obj(i, j) - (obj(i, k + 1) * cVl);
                            end
                        end
                    end
                end
            end
        end
        function b = solve(obj, b)
            iM = min(j + obj.lowerBw, obj.N - 1);
            for j = 0:(obj.N - 1)
                for i = (j + 1):iM
                    if obj(i, j + 1) ~= 0.0
                        b(i, :) = b(i, :) - obj(i, j + 1) * b(j + 1, :);
                    end
                end
            end
            for j = (obj.N - 1):-1:0
                b(j + 1, :) = b(j + 1, :) / obj(j + 1, j + 1);
                iM = max(0, j - obj.upperBw);
                for i = iM:(j - 1)
                    if obj(i + 1, j + 1) ~= 0.0
                        b(i + 1, :) = b(i + 1, :) - obj(i + 1, j + 1) * b(j + 1, :);
                    end
                end
            end
        end

        function b = solveAdj(obj, b)
            iM = min(j + obj.upperBw, obj.N - 1);
            for j = 0:(obj.N - 1)
                b(j + 1, :) = b(j + 1, :) / obj(j + 1, j + 1);
                for i = (j + 1):iM
                    if obj(j + 1, i + 1) ~= 0.0
                        b(i + 1, :) = b(i + 1, :) - obj(j + 1, i + 1) * b(j + 1, :);
                    end
                end
            end
            for j = (obj.N - 1):-1:0
                iM = max(0, j - obj.lowerBw);
                for i = iM:(j - 1)
                    if obj(j + 1, i + 1) ~= 0.0
                        b(i + 1, :) = b(i + 1, :) - obj(j + 1, i + 1) * b(j + 1, :);
                    end
                end
            end
        end
    end
end

