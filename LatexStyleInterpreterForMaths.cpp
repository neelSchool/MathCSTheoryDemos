#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cctype>
#include <algorithm>
#include <stdexcept>
#include <utility>
// Notes: There are a few cases that have not yet been adequately handled in this code (mainly where certain operations or values are undefined for example with division by zero). This code serves as a proof of concept 
// and is a continued work in progress which will be updated incrementally with additional edge case considerations (as was done thoroughly over time with the eccpointaddition script), and further support will be added
// for calculus (limits, differentiation and integration), and eventually complex number support and linear algebra operations (however these will take a long time to derive and implement efficiently). 
struct Expr;
using ExprPtr = std::shared_ptr<Expr>;

struct Expr {
    virtual ~Expr() = default;
    virtual std::string toLatex() const = 0;
    virtual ExprPtr simplify() const = 0;

    virtual bool isNumber(double& out) const {
        return false;
    }

    virtual std::string sortKey() const {
        return toLatex();
    }
};

ExprPtr num(double v);
ExprPtr sym(const std::string& name);
ExprPtr makeSum(std::vector<ExprPtr> terms);
ExprPtr makeProd(std::vector<ExprPtr> factors);
ExprPtr makePow(ExprPtr base, ExprPtr exp);
ExprPtr makeFunc(const std::string& name, ExprPtr arg);

bool isZero(const ExprPtr& e) {
    double v;
    return e->isNumber(v) && std::abs(v) < 1e-12;
}
bool isOne(const ExprPtr& e) {
    double v;
    return e->isNumber(v) && std::abs(v - 1.0) < 1e-12;
}
std::pair<double, ExprPtr> splitCoeff(const ExprPtr& e);
std::pair<ExprPtr, ExprPtr> splitPower(const ExprPtr& e);

struct Number : Expr {
    double val;
    explicit Number(double v)
        : val(std::abs(v) < 1e-12 ? 0.0 : v) {}
    bool isNumber(double& out) const override {
        out = val;
        return true;
    }
    std::string toLatex() const override {
        if (std::floor(val) == val) {
            return std::to_string(
                static_cast<long long>(val)
            );
        }
        std::string s = std::to_string(val);
        s.erase(
            s.find_last_not_of('0') + 1,
            std::string::npos
        );
        if (!s.empty() && s.back() == '.') {
            s.pop_back();
        }
        return s;
    }
    ExprPtr simplify() const override {
        return std::make_shared<Number>(val);
    }
    std::string sortKey() const override {
        return "0_num";
    }
};
struct Symbol : Expr {
    std::string name;
    explicit Symbol(std::string n)
        : name(std::move(n)) {}
    std::string toLatex() const override {
        return name;
    }
    ExprPtr simplify() const override {
        return std::make_shared<Symbol>(name);
    }
    std::string sortKey() const override {
        return "1_" + name;
    }
};
struct Sum;
struct Pow : Expr {
    ExprPtr base, exp;
    Pow(ExprPtr b, ExprPtr e)
        : base(std::move(b)), exp(std::move(e)) {}
    std::string toLatex() const override {
        std::string bStr = base->toLatex();
        if (bStr.size() > 1 && bStr[0] != '\\') {
            bStr = "(" + bStr + ")";
        }
        return bStr + "^{" + exp->toLatex() + "}";
    }
    ExprPtr simplify() const override {
        auto b = base->simplify();
        auto e = exp->simplify();
        double bv = 0, ev = 0;
        if (isZero(e)) return num(1);
        if (isOne(e)) return b;
        if (isZero(b)) return num(0);
        if (isOne(b)) return num(1);

        if (b->isNumber(bv) && e->isNumber(ev)) {
            return num(std::pow(bv, ev));
        }
        // Power raised to another power: (x^a)^b -> x^(a*b)
        if (auto innerPow =
            std::dynamic_pointer_cast<Pow>(b)) {
            return makePow(
                innerPow->base,
                makeProd({innerPow->exp, e})
            )->simplify();
        }
        return std::make_shared<Pow>(b, e);
    }
    std::string sortKey() const override {
        return base->sortKey() + "^" + exp->sortKey();
    }
};
struct Func : Expr {
    std::string name;
    ExprPtr arg;
    Func(std::string n, ExprPtr a)
        : name(std::move(n)), arg(std::move(a)) {}
    std::string toLatex() const override {
        return "\\" + name + "(" + arg->toLatex() + ")";
    }
    ExprPtr simplify() const override {
        auto a = arg->simplify();
        double av = 0;
        bool isNum = a->isNumber(av);
        std::string aTex = a->toLatex();
        // Natural Log identities
        if (name == "ln" || name == "log") {
            if (isOne(a)) return num(0);
            if (aTex == "e") return num(1);
            if (auto p = std::dynamic_pointer_cast<Pow>(a)) {
                if (p->base->toLatex() == "e") {
                    return p->exp;
                }
                return makeProd({
                    p->exp,
                    makeFunc(name, p->base)
                })->simplify();
            }
        }
        // Exponential identities
        if (name == "exp") {
            if (isZero(a)) return num(1);

            if (auto f = std::dynamic_pointer_cast<Func>(a)) {
                if (f->name == "ln") {
                    return f->arg;
                }
            }
        }
        // Trig exact vals
        if (name == "sin") {
            if (isZero(a) || aTex == "\\pi") {
                return num(0);
            }
        }
        if (name == "cos") {
            if (isZero(a)) return num(1);
            if (aTex == "\\pi") return num(-1);
        }
        if (name == "tan") {
            if (isZero(a) || aTex == "\\pi") {
                return num(0);
            }
        }
        if (isNum) {
            if (name == "sin") return num(std::sin(av));
            if (name == "cos") return num(std::cos(av));
            if (name == "tan") return num(std::tan(av));
            if (name == "ln" && av > 0) return num(std::log(av));
            if (name == "exp") return num(std::exp(av));
        }
        return std::make_shared<Func>(name, a);
    }
    std::string sortKey() const override {
        return "2_" + name + "(" + arg->sortKey() + ")";
    }
};
struct Product : Expr {
    std::vector<ExprPtr> factors;
    explicit Product(std::vector<ExprPtr> f)
        : factors(std::move(f)) {}
    std::string toLatex() const override {
        if (factors.empty()) return "1";
        std::string out;
        for (size_t i = 0; i < factors.size(); ++i) {
            double v = 0;
            bool isNum = factors[i]->isNumber(v);
            std::string tex = factors[i]->toLatex();
            if (i == 0 && isNum && v == -1.0 &&
                factors.size() > 1) {
                out += "-";
                continue;
            }
            if (std::dynamic_pointer_cast<Sum>(factors[i])) {
                tex = "(" + tex + ")";
            }
            if (i > 0 && !(i == 1 && out == "-")) {
                if (isNum || (!tex.empty() && tex[0] == '\\')) {
                    out += " \\cdot ";
                }
            }
            out += tex;
        }
        return out;
    }
    ExprPtr simplify() const override {
        std::vector<ExprPtr> flat;
        for (const auto& f : factors) {
            auto sf = f->simplify();
            if (auto subProd =
                std::dynamic_pointer_cast<Product>(sf)) {
                flat.insert(
                    flat.end(),
                    subProd->factors.begin(),
                    subProd->factors.end()
                );
            } else {
                flat.push_back(sf);
            }
        }
        double numCoeff = 1.0;
        std::map<
            std::string,
            std::pair<ExprPtr, std::vector<ExprPtr>>
        > groups;
        for (const auto& f : flat) {
            double v = 0;
            if (f->isNumber(v)) {
                if (std::abs(v) < 1e-12) return num(0);
                numCoeff *= v;
            } else {
                auto [base, exp] = splitPower(f);
                std::string key = base->toLatex();

                if (groups.find(key) == groups.end()) {
                    groups[key] = {base, {exp}};
                } else {
                    groups[key].second.push_back(exp);
                }
            }
        }
        std::vector<ExprPtr> result;
        if (std::abs(numCoeff - 1.0) > 1e-12) {
            result.push_back(num(numCoeff));
        }
        std::vector<ExprPtr> symTerms;
        for (auto& [key, pair] : groups) {
            ExprPtr totalExp =
                makeSum(pair.second)->simplify();

            if (!isZero(totalExp)) {
                symTerms.push_back(
                    makePow(pair.first, totalExp)->simplify()
                );
            }
        }
        std::sort(
            symTerms.begin(),
            symTerms.end(),
            [](const ExprPtr& a, const ExprPtr& b) {
                return a->sortKey() < b->sortKey();
            }
        );
        result.insert(
            result.end(),
            symTerms.begin(),
            symTerms.end()
        );
        if (result.empty()) return num(1);
        if (result.size() == 1) return result[0];
        return std::make_shared<Product>(result);
    }
};
struct Sum : Expr {
    std::vector<ExprPtr> terms;
    explicit Sum(std::vector<ExprPtr> t)
        : terms(std::move(t)) {}
    std::string toLatex() const override {
        if (terms.empty()) return "0";
        std::string out = terms[0]->toLatex();
        for (size_t i = 1; i < terms.size(); ++i) {
            std::string tStr = terms[i]->toLatex();

            if (!tStr.empty() && tStr[0] == '-') {
                out += " - " + tStr.substr(1);
            } else {
                out += " + " + tStr;
            }
        }
        return out;
    }
    ExprPtr simplify() const override {
        std::vector<ExprPtr> flat;
        for (const auto& t : terms) {
            auto st = t->simplify();
            if (auto subSum =
                std::dynamic_pointer_cast<Sum>(st)) {
                flat.insert(
                    flat.end(),
                    subSum->terms.begin(),
                    subSum->terms.end()
                );
            } else {
                flat.push_back(st);
            }
        }
        double constSum = 0.0;
        std::map<
            std::string,
            std::pair<ExprPtr, double>
        > likeTerms;
        for (const auto& t : flat) {
            double v = 0;

            if (t->isNumber(v)) {
                constSum += v;
            } else {
                auto [coeff, core] = splitCoeff(t);
                std::string key = core->toLatex();

                if (likeTerms.find(key) == likeTerms.end()) {
                    likeTerms[key] = {core, coeff};
                } else {
                    likeTerms[key].second += coeff;
                }
            }
        }

        std::vector<ExprPtr> result;

        for (auto& [key, pair] : likeTerms) {
            if (std::abs(pair.second) > 1e-12) {
                result.push_back(
                    makeProd({
                        num(pair.second),
                        pair.first
                    })->simplify()
                );
            }
        }
        std::sort(
            result.begin(),
            result.end(),
            [](const ExprPtr& a, const ExprPtr& b) {
                return a->sortKey() < b->sortKey();
            }
        );

        if (std::abs(constSum) > 1e-12) {
            result.push_back(num(constSum));
        }

        if (result.empty()) return num(0);
        if (result.size() == 1) return result[0];

        return std::make_shared<Sum>(result);
    }
};

ExprPtr num(double v) {
    return std::make_shared<Number>(v);
}
ExprPtr sym(const std::string& n) {
    return std::make_shared<Symbol>(n);
}
ExprPtr makeSum(std::vector<ExprPtr> t) {
    return std::make_shared<Sum>(std::move(t));
}
ExprPtr makeProd(std::vector<ExprPtr> f) {
    return std::make_shared<Product>(std::move(f));
}
ExprPtr makePow(ExprPtr b, ExprPtr e) {
    return std::make_shared<Pow>(
        std::move(b),
        std::move(e)
    );
}
ExprPtr makeFunc(const std::string& n, ExprPtr a) {
    return std::make_shared<Func>(
        n,
        std::move(a)
    );
}
std::pair<double, ExprPtr> splitCoeff(const ExprPtr& e) {
    if (auto prod = std::dynamic_pointer_cast<Product>(e)) {
        double c = 1.0;
        std::vector<ExprPtr> rest;

        for (const auto& f : prod->factors) {
            double v = 0;

            if (f->isNumber(v)) {
                c *= v;
            } else {
                rest.push_back(f);
            }
        }

        if (rest.empty()) return {c, num(1)};
        if (rest.size() == 1) return {c, rest[0]};

        return {c, makeProd(rest)};
    }

    return {1.0, e};
}
std::pair<ExprPtr, ExprPtr> splitPower(const ExprPtr& e) {
    if (auto p = std::dynamic_pointer_cast<Pow>(e)) {
        return {p->base, p->exp};
    }
    return {e, num(1)};
}
enum class TokType {
    NUM,
    ID,
    CMD,
    PLUS,
    MINUS,
    MUL,
    DIV,
    POW,
    LPAREN,
    RPAREN,
    LBRACE,
    RBRACE,
    END
};

struct Token {
    TokType type;
    std::string text;
};
class Lexer {
    std::string src;
    size_t pos = 0;
public:
    explicit Lexer(std::string s)
        : src(std::move(s)) {}
    std::vector<Token> tokenize() {
        std::vector<Token> tokens;
        while (pos < src.size()) {
            char c = src[pos];
            if (std::isspace(static_cast<unsigned char>(c))) {
                pos++;
                continue;
            }
            // Nums
            if (std::isdigit(static_cast<unsigned char>(c)) ||
                c == '.') {

                size_t start = pos;

                while (
                    pos < src.size() &&
                    (std::isdigit(static_cast<unsigned char>(src[pos])) ||
                     src[pos] == '.')
                ) {
                    pos++;
                }

                tokens.push_back({
                    TokType::NUM,
                    src.substr(start, pos - start)
                });
            }

            // LaTeX commands
            else if (c == '\\') {
                pos++;
                size_t start = pos;

                while (
                    pos < src.size() &&
                    std::isalpha(static_cast<unsigned char>(src[pos]))
                ) {
                    pos++;
                }

                std::string cmd =
                    src.substr(start, pos - start);

                if (cmd == "cdot") {
                    tokens.push_back({
                        TokType::MUL,
                        ""
                    });
                } else {
                    tokens.push_back({
                        TokType::CMD,
                        cmd
                    });
                }
            }
            // Vars
            else if (std::isalpha(static_cast<unsigned char>(c))) {
                tokens.push_back({
                    TokType::ID,
                    std::string(1, c)
                });

                pos++;
            }
            // Operators and brackets
            else {
                pos++;
                switch (c) {
                    case '+':
                        tokens.push_back({
                            TokType::PLUS, "+"
                        });
                        break;

                    case '-':
                        tokens.push_back({
                            TokType::MINUS, "-"
                        });
                        break;

                    case '*':
                        tokens.push_back({
                            TokType::MUL, "*"
                        });
                        break;

                    case '/':
                        tokens.push_back({
                            TokType::DIV, "/"
                        });
                        break;

                    case '^':
                        tokens.push_back({
                            TokType::POW, "^"
                        });
                        break;

                    case '(':
                        tokens.push_back({
                            TokType::LPAREN, "("
                        });
                        break;

                    case ')':
                        tokens.push_back({
                            TokType::RPAREN, ")"
                        });
                        break;

                    case '{':
                        tokens.push_back({
                            TokType::LBRACE, "{"
                        });
                        break;
                  
                    case '}':
                        tokens.push_back({
                            TokType::RBRACE, "}"
                        });
                        break;
                    default:
                        throw std::runtime_error(
                            std::string("Unexpected char: ") + c
                        );
                }
            }
        }
        tokens.push_back({
            TokType::END,
            ""
        });
        return tokens;
    }
};
class Parser {
    std::vector<Token> toks;
    size_t pos = 0;
    Token peek() const {
        return toks[pos];
    }
    Token consume() {
        return toks[pos++];
    }
    void expect(TokType t) {
        if (peek().type != t) {
            throw std::runtime_error(
                "Syntax error near: " + peek().text
            );
        }
        consume();
    }
    ExprPtr parseBraced() {
        expect(TokType::LBRACE);
        auto e = parseAddSub();
        expect(TokType::RBRACE);
        return e;
    }
public:
    explicit Parser(std::vector<Token> t)
        : toks(std::move(t)) {}
    // Add & subtract
    ExprPtr parseAddSub() {
        std::vector<ExprPtr> terms = {
            parseMulDiv()
        };
        while (
            peek().type == TokType::PLUS ||
            peek().type == TokType::MINUS
        ) {
            char op = consume().text[0];
            auto nextTerm = parseMulDiv();
            if (op == '-') {
                nextTerm = makeProd({
                    num(-1),
                    nextTerm
                });
            }
            terms.push_back(nextTerm);
        }
        return makeSum(terms);
    }
    // Mult & div
    ExprPtr parseMulDiv() {
        std::vector<ExprPtr> factors = {
            parsePower()
        };
        while (true) {
            if (peek().type == TokType::MUL) {
                consume();
                factors.push_back(parsePower());
            }
            else if (peek().type == TokType::DIV) {
                consume();
                factors.push_back(
                    makePow(parsePower(), num(-1))
                );
            }
            // Implicit multiplication:
            // 2x, 3sin(x), x(y+1)
            else if (
                peek().type == TokType::ID ||
                peek().type == TokType::LPAREN ||
                peek().type == TokType::CMD
            ) {
                factors.push_back(parsePower());
            }
            else {
                break;
            }
        }
        return makeProd(factors);
    }
    // Exponentiation
    ExprPtr parsePower() {
        auto base = parsePrimary();
        if (peek().type == TokType::POW) {
            consume();
            ExprPtr exp =
                (peek().type == TokType::LBRACE)
                ? parseBraced()
                : parsePrimary();
            return makePow(base, exp);
        }
        return base;
    }
    // Primary expr
    ExprPtr parsePrimary() {
        Token t = peek();
        // Unary minus
        if (t.type == TokType::MINUS) {
            consume();
            return makeProd({
                num(-1),
                parsePrimary()
            });
        }
        // Nums
        if (t.type == TokType::NUM) {
            consume();

            return num(std::stod(t.text));
        }
        // Vars
        if (t.type == TokType::ID) {
            consume();

            return sym(t.text);
        }
        // Parentheses
        if (t.type == TokType::LPAREN) {
            consume();
            auto e = parseAddSub();
            expect(TokType::RPAREN);
            return e;
        }
        // Braces
        if (t.type == TokType::LBRACE) {
            return parseBraced();
        }
        // LaTeX commands
        if (t.type == TokType::CMD) {
            consume();
            // Fractions: \frac{a}{b}
            if (t.text == "frac") {
                auto numerator = parseBraced();
                auto denominator = parseBraced();
                return makeProd({
                    numerator,
                    makePow(denominator, num(-1))
                });
            }

            if (t.text == "pi") {
                return sym("\\pi");
            }
            if (t.text == "sqrt") {
                return makePow(
                    parseBraced(),
                    num(0.5)
                );
            }
            // Functions: \sin(x), \cos(x), \ln(x), \exp(x)
            expect(TokType::LPAREN);
            auto arg = parseAddSub();
            expect(TokType::RPAREN);
            return makeFunc(t.text, arg);
        }
        throw std::runtime_error(
            "Unexpected token: " + t.text
        );
    }
};

// Simplify a LaTeX expression
std::string simplifyLatex(const std::string& input) {
    Lexer lex(input);
    Parser parser(lex.tokenize());

    return parser.parseAddSub()->simplify()->toLatex();
}
int main() {

    // Test expressions
    std::vector<std::string> tests = {
        "x + x",

        "3x + 5y - x + 2y + 4 - 1",

        "x \\cdot x^{2} \\cdot x^{3}",

        "\\frac{x^{5}}{x^{2}}",

        "(x^{2})^{3} + 2x^{6}",

        "\\ln(x^{3}) + 2\\ln(x)",

        "\\ln(e^{4x}) + \\exp(\\ln(2x))",

        "\\sin(0) + 5\\cos(0) + 3\\cos(\\pi)"
    };
    std::cout << "unit tests:\n\n";
    for (const auto& t : tests) {

        std::cout << "Input : " << t << "\n";

        try {
            std::cout << "Output: "
                      << simplifyLatex(t)
                      << "\n\n";
        }
        catch (const std::exception& e) {
            std::cout << "Error: "
                      << e.what()
                      << "\n\n";
        }
    }
    std::string line;
    std::cout << "Enter a LaTeX style expression to simplify.\n";
    std::cout << "Type 'quit' to quit.\n\n";
    std::cout << "> ";
    while (std::getline(std::cin, line) && line != "quit") {
        if (!line.empty()) {
            try {
                std::string result = simplifyLatex(line);

                std::cout << "= " << result << "\n";
            }
            catch (const std::exception& e) {
                std::cout << "Error: "
                          << e.what()
                          << "\n";
            }
        }
        std::cout << "> ";
    }
    std::cout << "Goodbye!\n";
    return 0;
}
/* 
compile and run command on my linux machine which uses g++
g++ -std=c++17 -Wall -Wextra -pedantic main.cpp -o simplifier
./simplifier
*/
// demo by Neel Mehendale
