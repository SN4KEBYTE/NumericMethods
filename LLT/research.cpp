#include "research.h"

void initial_testing()
{
    //for (int i = 1; i <= TESTS_NUM; i++)
    //{
    //    string end = to_string(i) + ".txt";
    //    ifstream in(INITIAL_TESTS_PATH + end);
    //    
    //    Matrix<type> m;
    //    vector<type> b;

    //    // input
    //    m.input(in, b);
    //    in.close();

    //    // solve
    //    m.factorization();
    //    m.forward(b);
    //    m.backward(b);

    //    // output
    //    ofstream out(INITIAL_RESULTS_PATH + end);
    //    out.precision(PREC);
    //    out.setf(std::ios::fixed);
    //    print_vector(b, out);
    //    out.close();

    //    // clear memory
    //    m.~Matrix();
    //    b.clear();
    //}
}

void ak_research()
{
    ofstream out(AK_RESULTS_PATH + "ak_res.csv");
    out.setf(std::ios::fixed);
    out.precision(numeric_limits<double>::max_digits10);
    out.imbue(locale(""));

    ifstream in(AK_TESTS_PATH + "ak.txt");

    auto es = get_exact_solution(AK_DIM);

    for (int i = 0; i < AK_NUM; i++)
    {
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

        float coef_f = pow(10, -i);
        m_float.ak_add(coef_f);

        m_float.factorization();
        m_float.forward(b_float);
        m_float.backward(b_float);

        // solve for double using regular factorization

        double coef_d = pow(10, -i);
        m_double.ak_add(coef_d);

        m_double.factorization();
        m_double.forward(b_double);
        m_double.backward(b_double);

        // solve for float using factorization_d

        m_float_d.ak_add(coef_f);

        m_float_d.factorization_d();
        m_float_d.forward(b_float_d);
        m_float_d.backward(b_float_d);

        dump_research(out, i, b_float, b_float_d, b_double, es);
    }
}

void gilbert_research()
{
    // TODO
    //for (int i = 0; i < GILBERT_NUM; i++)
    //{
    //    string end = "/matrixG" + to_string(i + 2) + ".txt";

    //    ifstream in(GILBERT_TESTS_PATH + end);
    //    Matrix<type> m;
    //    vector<type> b;

    //    // input
    //    m.input(in, b);
    //    in.close();

    //    // solve
    //    m.factorization();
    //    m.forward(b);
    //    m.backward(b);

    //    // b = xk
    //    auto x_ast = get_x_asterisk(m.getDimension());

    //    cout << "K = " << i + 2 << endl;

    //    cout << "xk = ";
    //    ShowVector(b, cout);

    //    auto diff = vec_diff(x_ast, b);
    //    cout << "diff = ";
    //    ShowVector(diff, cout);

    //    cout << "(xk, xk) = " << scalar_product(b, b) << endl;
    //    cout << "(x*-xk, x*-xk) = " << scalar_product(diff, diff) << endl;
    //    cout << "----------------------------------------" << endl;

    //    // output and clear memory.
    //    ofstream out(GILBERT_RESULTS_PATH + TYPE_STR + end);
    //    out.precision(PREC);
    //    out.setf(std::ios::fixed);
    //    ShowVector<type>(b, out);
    //    out.close();

    //    m.~Matrix();
    //    b.clear();
    //}
}

vector<int> get_exact_solution(const int &N)
{
    vector<int> x(N);

    for (size_t i = 0; i < x.size(); i++)
        x[i] = i + 1;

    return x;
}