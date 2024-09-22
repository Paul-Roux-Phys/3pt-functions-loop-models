#include "arsnsym.h"
#include "arscomp.h"

template <typename W, typename VP>
class Matrix {
private:
    VP* _v;
    int _size;
    // Key* _state_table;
    int _nev;
    int _nconv;
    #if (MATRIX_TYPE == REAL)
    ARNonSymStdEig<W, Matrix> _eig_prob;
    #else
    ARCompStdEig<W, Matrix> _eig_prob;
    #endif
    Real *_eig_r, *_eig_i;

public:
    Key* _state_table;

    Matrix(VP* v, int nev=2, const char* which="LM") : _nev(nev) {
        v->generate_space();
        _v = v;
        _size = v->size();
        make_state_table();
        _v->clear();
        _eig_prob.DefineParameters(size(), nev, this, &Matrix::MultMx, which);
    }

    ~Matrix() {
        delete[] _state_table;
    }

    int size() { return _size; }

    void MultMx(W* x, W* y) {
        // place x in the hash table
        for (int i=0; i < _size; i++) {
            (*_v->getempty())[_state_table[i]] = x[i];
        }
        _v->swap_data_position();

        // multiply by the transfer matrix
        _v->transfer();

        // put the result back into y
        for (int i=0; i < _size; i++) {
            y[i] = (*_v)[_state_table[i]];
        }

        _v->clear();
    }

    void make_state_table() {
        _state_table = new Key[_size];
        int i=0;
        for (auto it = _v->get()->draining_begin(); it != _v->get()->draining_end(); it++) {
            BasisVector<W> b(*it);
            _state_table[i] = b.key;
            i++;
        }
    }

    void print_state_table() {
        for (int i=0; i < _size; i++) {
            cout << std::setw(3) << i+1 << ": " << _state_table[i] << endl;
        }
    }

    W operator()(int i, int j) {
        double* x = new W[_size];
        double* y = new W[_size];

        for (int k=0; k < _size; k++) { x[k] = 0; }
        x[i] = 1;

        MultMx(x, y);

        W res = y[j];
        
        delete[] x;
        delete[] y;
        return res;
    }

    void print(int column_width = 10) {
        cout << "Matrix of dimension " << _size << ":" << endl;
        cout << "[" << endl;;
        for (int i=0; i < _size; i++) {
            cout << std::setw(4) << " ";
            for (int j=0; j < _size; j++) {
                cout << std::left << std::setw(column_width) << (*this)(i, j);
            }
            if (i != _size-1) {
                cout << ";" << endl;
            }
        }
        cout << endl << "]" << endl;
    }

    bool isapprox(Real a, Real b, Real eps) {
        return (abs(a - b) < eps || abs(1-a/b) < eps);
    }

    void remove_duplicate_eigvals(Real eps=1e-10) {
        Real* tmp_r = new Real[_nconv];
        Real* tmp_i = new Real[_nconv];
        int j=0;
        for (int i = 0; i < _nconv-1; i++) {
            if (!isapprox(_eig_r[i], _eig_r[i+1], eps) \
                || !isapprox(_eig_i[i], _eig_i[i+1], i+1)) 
            {
                tmp_r[j] = _eig_r[i];
                tmp_i[j] = _eig_i[i];
                j++;
            }
        }
        tmp_r[j] = _eig_r[_nconv-1];
        tmp_i[j] = _eig_i[_nconv-1];
        j++;

        delete[] _eig_r;
        delete[] _eig_i;
        _nconv = j;
        _eig_r = new Real[_nconv];
        _eig_i = new Real[_nconv];

        for (int i = 0; i < j; i++) {
            _eig_r[i] = tmp_r[i];
            _eig_i[i] = tmp_i[i];
        }
        
        delete[] tmp_r;
        delete[] tmp_i;
    }

    void find_eigenvalues(bool no_mult=false, Real eps = 1e-10) {
        _nconv = _eig_prob.FindEigenvalues();
        _eig_r = new double[_nconv];
        _eig_i = new double[_nconv];
        _eig_prob.Eigenvalues(_eig_r, _eig_i);
        if (no_mult)
        {
            remove_duplicate_eigvals(eps);
        }

    }

    std::pair<double*, double*> get_eigenvalues() {
        return std::make_pair(_eig_r, _eig_i);
    }

    // Print the eigenvalues to a string. If nev is set to zero, then all computed eigenvalues are printed.
    std::string sprint_eigenvalues(int nev=0, int width_real_part=20, int width_imag_part=10) {
        std::ostringstream oss;

        if (nev == 0) {
            nev = _nconv;
        }

        if (nev > _nconv) {
            nev = _nconv;
        }

        oss << std::setprecision(std::numeric_limits<Real>::max_digits10);
        for (int i = _nconv-1; i >= _nconv-nev; i--) {
            oss << "(" << std::left << std::setw(width_real_part) << _eig_r[i] << ", "
                             << std::setw(width_imag_part) << _eig_i[i] << ")" << endl;
        }

        return oss.str();
    }
};