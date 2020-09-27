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
    //ofstream out(AK_RESULTS_PATH + TYPE_STR + ".csv");
    //out.precision(PREC);
    //out.setf(std::ios::fixed);
    //out.imbue(locale(""));

    //for (int i = 0; i < AK_NUM; i++)
    //{
    //    ifstream in(AK_TESTS_PATH + "/matrixAK" + to_string(i) + ".txt");
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
    //    // TODO

    //    m.~Matrix();
    //    b.clear();
    //}
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