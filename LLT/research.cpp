#include "research.h"

void initial_testing()
{
    for (int i = 1; i <= TESTS_NUM; i++)
    {
        string end = to_string(i) + ".txt";
        ifstream in(INITIAL_TESTS_PATH + end);
        
        Matrix<type> m;
        vector<type> b;

        // input
        m.input(in, b);
        in.close();

        // solve
        m.factorization(m);
        m.forward(b, b);
        m.backward(b, b);

        // output
        ofstream out(INITIAL_RESULTS_PATH + end);
        out.precision(PREC);
        out.setf(std::ios::fixed);
        print_vector(b, out);
        out.close();

        // clear memory
        m.~Matrix();
        b.clear();
    }
}

void ak_research()
{
    ofstream out(AK_RESULTS_PATH + "ak_res.csv");
    out.setf(std::ios::fixed);
    out.precision(numeric_limits<double>::max_digits10);
    out.imbue(locale(""));

    ifstream in(AK_TESTS_PATH + "ak.txt");

    auto es = get_exact_solution<int>(AK_DIM);

    for (int i = 0; i < AK_NUM; i++)
    {
        Matrix<float> m_float;
        vector<float> b_float;
        
        Matrix<float> m_float_d;
        vector<float> b_float_d;

        Matrix<double> m_double;
        vector<double> b_double;

        m_float.input(in, b_float);

        in.seekg(0, ios::beg);
        m_double.input(in, b_double);

        in.seekg(0, ios::beg);
        m_float_d.input(in, b_float_d);

        in.seekg(0, ios::beg);

        // solve for float using regular factorization
        float coef_f = pow(10, -i);
        
        m_float.ak_add(coef_f);
        b_float[0] += coef_f;

        m_float.factorization(m_float);
        m_float.forward(b_float, b_float);
        m_float.backward(b_float, b_float);

        // solve for double using regular factorization
        double coef_d = pow(10, -i);
        
        m_double.ak_add(coef_d);
        b_double[0] += coef_d;

        m_double.factorization(m_double);
        m_double.forward(b_double, b_double);
        m_double.backward(b_double, b_double);

        // solve for float using factorization_d
        m_float_d.ak_add(coef_d);
        b_float_d[0] += coef_d;

        m_float_d.factorization_d(m_float_d);
        m_float_d.forward_d(b_float_d, b_float_d);
        m_float_d.backward(b_float_d, b_float_d);

        dump_research(out, i, b_float, b_float_d, b_double, es);
    }
}

void gilbert_research()
{
   ofstream out(GILBERT_RESULTS_PATH + "gilbert_res.csv");
   out.setf(std::ios::fixed);
   out.precision(numeric_limits<double>::max_digits10);
   out.imbue(locale(""));

   for (size_t i = 2; i <= GILBERT_NUM; i++)
   {
      ifstream in(GILBERT_TESTS_PATH + "gilbert" + to_string(i) + ".txt");
      auto es = get_exact_solution<int>(i);

      Matrix<float> m_float;
      Matrix<float> m_float_d;
      vector<float> b_float;
      vector<float> b_float_d;

      Matrix<double> m_double;
      vector<double> b_double;

      m_float.input(in, b_float);

      in.seekg(0, ios::beg);
      m_double.input(in, b_double);

      in.seekg(0, ios::beg);
      m_float_d.input(in, b_float_d);

      in.seekg(0, ios::beg);

      // solve for float using regular factorization
      m_float.factorization(m_float);
      m_float.forward(b_float, b_float);
      m_float.backward(b_float, b_float);

      // solve for double using regular factorization
      m_double.factorization(m_double);
      m_double.forward(b_double, b_double);
      m_double.backward(b_double, b_double);

      // solve for float using factorization_d
      m_float_d.factorization_d(m_float_d);
      m_float_d.forward(b_float_d, b_float_d);
      m_float_d.backward(b_float_d, b_float_d);

      dump_research(out, i, b_float, b_float_d, b_double, es);
   }
}

void gaussian_research()
{
    vector<vector<double>> md(G_DIM);
    vector<double> b(G_DIM);

    for (auto &row : md)
        row = vector<double>(G_DIM + 1);

    Matrix<double> mp;

    ifstream in(GAUSSIAN_TESTS_PATH + "dense.txt");

    for (size_t i = 0; i < md.size(); i++)
        for (size_t j = 0; j < md[0].size(); j++)
            in >> md[i][j];

    in.close();
    in.open(GAUSSIAN_TESTS_PATH + "profile.txt");

    mp.input(in, b);
    auto b_copy = b;
    gauss(md, b_copy);
    
    mp.factorization(mp);
    mp.forward(b, b);
    mp.backward(b, b);

    auto es = get_exact_solution<double>(AK_DIM);

    ofstream out(GAUSSIAN_RESULTS_PATH + "gaussian_res.csv");
    out.setf(std::ios::fixed);
    out.precision(numeric_limits<double>::max_digits10);
    out.imbue(locale(""));
    dump_gaussian_research(out, es, b, b_copy);
}

template <typename T>
int gauss(std::vector<std::vector<T>> &A, std::vector<T> &x)
{
    int n = (int)A.size();
    int m = (int)A[0].size() - 1;

    std::vector<int> where(m, -1);

    for (int col = 0, row = 0; col < m && row < n; ++col)
    {
        int sel = row;

        for (int i = row; i < n; ++i)
            if (abs(A[i][col]) > abs(A[sel][col]))
                sel = i;

        if (abs(A[sel][col]) < EPS)
            continue;

        for (int i = col; i <= m; ++i)
            swap(A[sel][i], A[row][i]);

        where[col] = row;

        for (int i = 0; i < n; ++i)
            if (i != row)
            {
                double c = A[i][col] / A[row][col];

                for (int j = col; j <= m; ++j)
                    A[i][j] -= A[row][j] * c;
            }

        ++row;
    }

    x.assign(m, 0);

    for (int i = 0; i < m; ++i)
        if (where[i] != -1)
            x[i] = A[where[i]][m] / A[where[i]][i];

    for (int i = 0; i < n; ++i)
    {
        double sum = 0;

        for (int j = 0; j < m; ++j)
            sum += x[j] * A[i][j];

        if (abs(sum - A[i][m]) > EPS)
            return 0;
    }

    for (int i = 0; i < m; ++i)
        if (where[i] == -1)
            return INFINITY;

    return 1;
}

template<typename T>
vector<T> get_exact_solution(const int &N)
{
    vector<T> x(N);

    for (size_t i = 0; i < x.size(); i++)
        x[i] = i + 1;

    return x;
}

template<typename T>
void print_vector(std::vector<T> &v, std::ostream &out)
{
    for (size_t i = 0; i < v.size(); i++)
        out << v[i] << std::endl;
}