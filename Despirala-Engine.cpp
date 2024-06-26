#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <random>
#include <math.h>
#include <time.h>

const double EPS = 1e-6;
const double INF = 1e6;

// Game configuration
const int NUM_DICE = 6;
const int NUM_SIDES = 6;
const int NUM_COMBOS = 14;
const int GOODS_PER_TURN = 5;
const int POINTS_PER_GOOD = 1;

// Depends on game configuration
const int MAX_ROLL = 46656; // NUM_DICE ^ NUM_SIDES
const int MAX_CODE = 1 << (NUM_SIDES + NUM_DICE);
const int DIFF_CODES = 30; // Sum_{0 <= k <= NUM_DICE} #ways to partition k into at most NUM_SIDES partitions
const int NUM_MASKS = 1 << NUM_COMBOS;
const int MAX_GOODS = NUM_COMBOS * GOODS_PER_TURN;

const int NUM_TRIALS_1 = 100000; // Trials for rolls distribution
int NUM_TRIALS_2; // Trials per state

typedef std::array<int, NUM_SIDES> Occurs;

struct Target
{
    const std::vector<int> args;
    const int points;
    const Occurs occurs;

    Target(int points, const Occurs& occurs):
        args({}),
        points(points),
        occurs(occurs) {}
    Target(const std::vector<int>& args, int points, const Occurs& occurs):
        args(args),
        points(points),
        occurs(occurs) {}
};

struct Combo
{
    Combo(const std::string& name, int points):
        name(name),
        points(points) {}
    std::string getName() const
    {
        return name;
    }
    virtual int getNumArgs() const
    {
        return 0;
    }
    virtual int getCollectNumber() const
    {
        return 0;
    }
    virtual std::vector<Target> getTargets(const Occurs& diceOccurs)
    {
        return {};
    }
    virtual Target getTargetByArgs(const std::vector<int>& args) const
    {
        return Target(0, {});
    }

protected:
    const std::string name;
    const int points;

    static Occurs sortedOccurs(Occurs occurs)
    {
        std::sort(occurs.begin(), occurs.end());
        return occurs;
    }
};

struct CollectCombo : Combo
{
    CollectCombo(const std::string& name, int number):
        Combo(name, 0),
        number(number) {}
    int getNumArgs() const
    {
        return 1;
    }
    int getCollectNumber() const
    {
        return number;
    }
    

protected:
    int number;
};

struct FixedCombo : Combo
{
public:
    FixedCombo(const std::string& name, int points, const Occurs& occurs):
        Combo(name, points),
        occurs(occurs) {}
    std::vector<Target> getTargets(const Occurs& diceOccurs)
    {
        return {Target(points, occurs)};
    }
    Target getTargetByArgs(const std::vector<int>& args) const
    {
        return Target(points, occurs);
    }

protected:
    const Occurs occurs;
};

struct PermCombo : Combo
{
public:
    PermCombo(const std::string& name, int points, const Occurs& occurs):
        Combo(name, points),
        templateOccurs(sortedOccurs(occurs)) {}
    int getNumArgs() const
    {
        int numArgs = 0;
        for (int occ : templateOccurs)
        {
            if (occ > 0) ++numArgs;
        }
        return numArgs;
    }
    std::vector<Target> getTargets(const Occurs& diceOccurs)
    {
        int code = 0;
        for (int diceOccur : diceOccurs)
        {
            code = code * (NUM_DICE + 1) + diceOccur;
        }
        if (cache.find(code) == cache.end())
        {
            cache[code] = rawGetTargets(diceOccurs);
        }
        return cache[code];
    }
    Target getTargetByArgs(const std::vector<int>& args) const
    {
        Occurs occurs;
        makeOccursByArgs(templateOccurs, args, occurs);

        return Target(args, points, occurs);
    }

protected:
    const Occurs templateOccurs;
    std::unordered_map<int, std::vector<Target>> cache;

    virtual std::vector<Target> rawGetTargets(const Occurs& diceOccurs) const
    {
        std::array<std::pair<int, int>, NUM_SIDES> diceOccursPairs;
        makePairs(diceOccurs, diceOccursPairs);
        std::sort(diceOccursPairs.begin(), diceOccursPairs.end());

        Occurs occurs;
        std::vector<int> args = makeOccurs(templateOccurs, diceOccursPairs, occurs);

        return {Target(args, points, occurs)};
    }
    static void makePairs(const Occurs& diceOccurs, std::array<std::pair<int, int>, NUM_SIDES>& diceOccursPairs)
    {
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            diceOccursPairs[i] = {diceOccurs[i], i};
        }
    }
    static std::vector<int> makeOccurs(const Occurs& templateOccurs, const std::array<std::pair<int, int>, NUM_SIDES>& diceOccursPairs, Occurs& occurs)
    {
        std::vector<int> args;
        for (int i = NUM_SIDES - 1; i >= 0; --i)
        {
            occurs[diceOccursPairs[i].second] = templateOccurs[i];
            if (templateOccurs[i] > 0) args.push_back(diceOccursPairs[i].second + 1);
        }
        return args;
    }
    static void makeOccursByArgs(const Occurs& templateOccurs, const std::vector<int>& args, Occurs& occurs)
    {
        occurs.fill(0);
        for (int i = 0; i < (int) args.size(); ++i)
        {
            occurs[args[i] - 1] = templateOccurs[NUM_SIDES - i - 1];
        }
    }
};

struct SPermCombo : PermCombo
{
    SPermCombo(const std::string& name, const Occurs& occurs):
        PermCombo(name, 0, occurs) {};
    Target getTargetByArgs(const std::vector<int>& args) const
    {
        Occurs occurs;
        makeOccursByArgs(templateOccurs, args, occurs);
        int currPoint = evalOccurs(occurs);

        return Target(args, currPoint, occurs);
    }

protected:
    std::vector<Target> rawGetTargets(const Occurs& diceOccurs) const
    {
        std::array<std::pair<int, int>, NUM_SIDES> diceOccursPairs;
        makePairs(diceOccurs, diceOccursPairs);
        std::sort(diceOccursPairs.begin(), diceOccursPairs.end());

        Occurs currTemplate(templateOccurs);

        int maxPoints = -1;
        std::vector<Target> ts;
        Occurs occurs;
        do
        {
            std::vector<int> args = makeOccurs(currTemplate, diceOccursPairs, occurs);
            int currPoints = evalOccurs(occurs);
            if (currPoints > maxPoints)
            {
                maxPoints = currPoints;
                ts.push_back(Target(args, currPoints, occurs));
            }
        }
        while (std::next_permutation(currTemplate.begin(), currTemplate.end()));

        return ts;
    }
    static int evalOccurs(const Occurs& occurs)
    {
        int points = 0;
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            points += (i + 1) * occurs[i];
        }
        return points;
    }
};

// Game configuration
Combo* combos[]={new CollectCombo("Collect", 1),
                 new CollectCombo("Collect", 2),
                 new CollectCombo("Collect", 3),
                 new CollectCombo("Collect", 4),
                 new CollectCombo("Collect", 5),
                 new CollectCombo("Collect", 6),
                 new   SPermCombo("Three pairs",        {2, 2, 2, 0, 0, 0}),
                 new   SPermCombo("Two triples",        {3, 3, 0, 0, 0, 0}),
                 new    PermCombo("Four of a kind", 40, {4, 0, 0, 0, 0, 0}),
                 new   FixedCombo("Kamerun",        45, {0, 0, 0, 1, 2, 3}),
                 new   FixedCombo("Straight",       50, {1, 1, 1, 1, 1, 1}),
                 new    PermCombo("Six of a kind",  60, {6, 0, 0, 0, 0, 0}),
                 new   FixedCombo("General",        70, {0, 0, 0, 0, 0, 6}),
                 new   FixedCombo("Despirala",      80, {5, 0, 0, 0, 0, 1})};

// Depends on game configuration
const int MIN_COMB_DICE = 4;

int occursToCode(Occurs occurs)
{
    std::sort(occurs.begin(), occurs.end(), std::greater<int>());
    int code = 1;
    int curr = 0;
    while (curr < NUM_SIDES)
    {
        while (occurs[curr])
        {
            --occurs[curr];
            code *= 2;
        }
        code = code * 2 + 1;
        ++curr;
    }
    code /= 2;
    return code;
}

bool codeToOccurs(int code, Occurs& occurs)
{
    occurs.fill(0);
    int n = NUM_SIDES + NUM_DICE;
    int instrNum = 0;
    int instrs[n + 1];
    while (code && instrNum <= n)
    {
        instrs[instrNum] = code % 2;
        code /= 2;
        ++instrNum;
    }
    if (instrNum == 0 || instrNum > n) return false;
    --instrNum;
    int curr = 0;
    while (instrNum && curr < NUM_SIDES)
    {
        --instrNum;
        if (instrs[instrNum]) ++curr;
        else ++occurs[curr];
    }
    if (curr != NUM_SIDES - 1) return false;
    for (int i = 1; i < NUM_SIDES; ++i)
    {
        if (occurs[i - 1] < occurs[i]) return false;
    }
    return true;
}

void rollToOccurs(int roll, Occurs& occurs)
{
    occurs.fill(0);
    while (roll)
    {
        ++occurs[roll % NUM_SIDES];
        roll /= NUM_SIDES;
    }
}

std::mt19937 generator;
std::uniform_int_distribution<int> diceDistribution(0, NUM_SIDES - 1);

int randDiceRoll()
{
    return diceDistribution(generator);
}

void randOcc(Occurs& occurs)
{
    occurs.fill(0);
    for (int i = 0; i < NUM_DICE; ++i)
    {
        ++occurs[randDiceRoll()];
    }
}

bool randRemOcc(int extraDice, Occurs& occurs)
{
    int numFreeDice = std::accumulate(occurs.begin(), occurs.end(), extraDice);
    bool changed = false;
    for (int i = 0; i < numFreeDice; ++i)
    {
        int side = randDiceRoll();
        if (occurs[side])
        {
            --occurs[side];
            changed = true;
        }
    }
    return changed;
}

void remOcc(const Occurs& diceOccurs, Occurs& occurs)
{
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        occurs[i] = std::max(occurs[i] - diceOccurs[i], 0);
    }
}

std::unordered_map<int, int> codeToIndex;
int indexToCode[DIFF_CODES];

void generateLefts()
{
    Occurs occurs;
    int lastIdx = 0;
    for (int i = 0; i <= MAX_CODE; ++i)
    {
        if (codeToOccurs(i, occurs))
        {
            indexToCode[lastIdx] = i;
            codeToIndex[i] = lastIdx;
            ++lastIdx;
        }
    }
}

const int MAX_EXTRA_DICE = NUM_DICE - MIN_COMB_DICE;
double rollsDistr[MAX_EXTRA_DICE + 1][DIFF_CODES][MAX_GOODS + 1];

void findRollsDistrSingle(int ed, int idx)
{
    Occurs occurs;
    Occurs occurs2;
    int code = indexToCode[idx];
    codeToOccurs(code, occurs);
    for (int t = 0; t < NUM_TRIALS_1; ++t)
    {
        occurs2 = occurs;
        bool changed = false;
        int rolls = 0;
        while (!changed)
        {
            changed = randRemOcc(ed, occurs2);
            ++rolls;
        }
        int newCode = occursToCode(occurs2);
        int newIdx = codeToIndex[newCode];
        for (int i = 0; i <= MAX_GOODS - rolls; ++i)
        {
            rollsDistr[ed][idx][i + rolls] += rollsDistr[ed][newIdx][i] / NUM_TRIALS_1;
        }
    }
}

void findRollsDistr()
{
    for (int ed = 0; ed <= MAX_EXTRA_DICE; ++ed)
    {
        rollsDistr[ed][0][0] = 1;
        for (int i = 1; i < DIFF_CODES; ++i)
        {
            findRollsDistrSingle(ed, i);
        }
    }
}

double leftDistr[NUM_DICE + 1][NUM_DICE + 1];

int fact(int n)
{
    return n ? n * fact(n - 1) : 1;
}

void findLeftDistr()
{
    for (int i = 0; i <= NUM_DICE; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            leftDistr[i][j] = pow(5, j) / pow(6, i) * fact(i) / (fact(j) * fact(i - j));
        }
    }
}

#define M_STOP_COLL -1
#define M_CONT_COLL -2
#define M_REROLL -3
#define M_ERROR -100

struct Move
{
    int id;
    std::vector<int> args;
    double score;

    Move(int id, const std::vector<int>& args, double score):
        id(id),
        args(args),
        score(score) {}
    
    Move(int id, double score):
        Move(id, {}, score) {}
    
    Move(const std::string& name, const std::vector<int>& args);
    
    std::string toString() const;
};

bool operator<(const Move& a, const Move& b)
{
    return a.score < b.score;
}

bool isFree(int free, int idx)
{
    return (free >> idx) & 1;
}

int setUsed(int free, int idx)
{
    return free - (1 << idx);
}

bool isFound[NUM_MASKS][MAX_GOODS + 1];
double score[NUM_MASKS][MAX_GOODS + 1];
bool isFoundColl[NUM_MASKS][MAX_GOODS + 1][NUM_SIDES][NUM_DICE + 1];
double collScore[NUM_MASKS][MAX_GOODS + 1][NUM_SIDES][NUM_DICE + 1];

double getScore(int free, int goods);
double getScoreCont(int free, int goods);
double getScoreColl(int free, int goods, int num, int left);

double getScoreCollContinue(int free, int goods, int num, int left)
{
    double sc = 0;
    for (int i = 0 ; i <= left; ++i)
    {
        sc += leftDistr[left][i] * ((left - i) * num + getScoreColl(free, goods - 1, num, i));
    }
    return sc;
}

Move getMoveColl(int free, int goods, int num, int left)
{
    Move best = Move(M_STOP_COLL, getScore(free, goods));
    if (goods && left)
    {
        Move option = Move(M_CONT_COLL, getScoreCollContinue(free, goods, num, left));
        best = std::max(best, option);
    }
    return best;
}

double getScoreColl(int free, int goods, int num, int left)
{
    if (!isFoundColl[free][goods][num][left])
    {
        isFoundColl[free][goods][num][left] = true;
        collScore[free][goods][num][left] = getMoveColl(free, goods, num, left).score;
    }
    return collScore[free][goods][num][left];
}

int calcExtraDice(const Occurs& occurs)
{
    return NUM_DICE - std::accumulate(occurs.begin(), occurs.end(), 0);
}

double getMoveScore(int free, int goods, const Occurs& occurs, int extraDice, int reward)
{
    double score = 0;
    int code = occursToCode(occurs);
    int idx = codeToIndex[code];

    double succP = 0;
    for (int i = 0; i <= goods; ++i)
    {
        double currP = rollsDistr[extraDice][idx][i];
        if (currP > 0)
        {
            succP += currP;
            score += currP * getScore(free, goods - i);
        }
    }
    score += succP * reward;
    score += (1 - succP) * getScore(free, 0);
    return score;
}

Move getMove(int free, int goods, const Occurs& diceOccurs)
{
    Move best = Move(M_ERROR, -INF);
    for (int i = 0; i < NUM_COMBOS; ++i)
    {
        if (!isFree(free, i)) continue;
        int newFree = setUsed(free, i);
        int num = combos[i]->getCollectNumber();
        if (num)
        {
            int curr = diceOccurs[num - 1];
            Move option = Move(i, {num}, num * curr + getScoreColl(newFree, goods, num, NUM_DICE - curr));
            best = std::max(best, option);
        }
        else
        {
            std::vector<Target> ts = combos[i]->getTargets(diceOccurs);
            for (const Target& t : ts)
            {
                Occurs newOccurs = t.occurs;
                int extraDice = calcExtraDice(newOccurs);
                remOcc(diceOccurs, newOccurs);
                Move option = Move(i, t.args, getMoveScore(newFree, goods, newOccurs, extraDice, t.points));
                best = std::max(best, option);
            }
        }
    }
    if (goods > 0)
    {
        Move option = Move(M_REROLL, getScoreCont(free, goods - 1));
        best = std::max(best, option);
    }
    return best;
}

int statesVisCnt = 0;

double getScoreCont(int free, int goods)
{
    if (isFound[free][goods]) return score[free][goods];
    isFound[free][goods] = true;
    if (free == 0) score[free][goods] = goods * POINTS_PER_GOOD;
    else
    {
        Occurs diceOccurs;
        for (int t = 0; t < NUM_TRIALS_2; ++t)
        {
            randOcc(diceOccurs);
            Move mov = getMove(free, goods, diceOccurs);
            score[free][goods] += mov.score / NUM_TRIALS_2;
        }
    }
    ++statesVisCnt;
    if (statesVisCnt % 10000 == 0) std::cerr << ".";
    return score[free][goods];
}

double getScore(int free, int goods)
{
    return free ? getScoreCont(free, goods + GOODS_PER_TURN) : getScoreCont(free, goods);
}

double getInitialScore()
{
    return getScore(NUM_MASKS - 1, 0);
}

std::string stringToLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](char c){ return std::tolower(c); });
    return s;
}

bool isNumber(const std::string& s)
{
    return all_of(s.begin(), s.end(), [](char c){ return std::isdigit(c); });
}

std::string Move::toString() const
{
    std::string s;
    switch (id)
    {
    case M_ERROR:
        s = "Error";
        break;

    case M_STOP_COLL:
        s = "Stop collecting";
        break;

    case M_CONT_COLL:
        s = "Continue collecting";
        break;

    case M_REROLL:
        s = "Reroll";
        break;

    default:
        s = combos[id]->getName();
    }
    for (int a : args)
    {
        s += " " + std::to_string(a);
    }
    return s;
}

Move::Move(const std::string& name, const std::vector<int>& args):
    Move(M_ERROR, args, 0)
{
    std::unordered_set<int> freeArgs;
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        freeArgs.insert(i + 1);
    }

    for (int arg : args)
    {
        auto it = freeArgs.find(arg);
        if (it == freeArgs.end()) return;
        freeArgs.erase(it);
    }

    if (name == "stop collecting" || name == "stop") id = M_STOP_COLL;
    else if (name == "continue collecting" || name == "continue") id = M_CONT_COLL;
    else if (name == "reroll") id = M_REROLL;
    else
    {
        for (int i = 0; i < NUM_COMBOS; ++i)
        {
            if (name == stringToLower(combos[i]->getName()) && (int) args.size() == combos[i]->getNumArgs())
            {
                int num = combos[i]->getCollectNumber();
                if (!num || (!args.empty() && num == args[0]))
                {
                    id = i;
                    break;
                }
            }
        }
    }
}

const int NUM_REASONS = 3;

#define R_LUCK 0
#define R_MISTAKE 1
#define R_NONE 2

struct State
{
    int free;
    int goods;
    int score;
    double expScore;
    int fail;

    double scoreByReason[NUM_REASONS];

    State():
        free(NUM_MASKS - 1),
        goods(0),
        score(0),
        expScore(getInitialScore()),
        fail(false)
    {
        std::fill(scoreByReason, scoreByReason + NUM_REASONS, 0);
    }

    int getTurn()
    {
        int turn = 1;
        for (int i = 0; i < NUM_COMBOS; ++i)
        {
            if ((free >> i & 1) == 0) ++turn;
        }
        return turn;
    }

    void updateExpScore(double newExpScore, int reason, bool printEvalLM, const std::string& bestMvName="")
    {
        newExpScore += score;
        double delta = newExpScore - expScore;
        bool print = printEvalLM && (reason == R_LUCK || fabs(newExpScore - expScore) > EPS);
        if (print)
        {
            if (reason == R_LUCK && delta >= 0) std::cout << "Good luck: ";
            else if (reason == R_LUCK && delta < 0) std::cout << "Bad luck: ";
            else if (reason == R_MISTAKE && fabs(delta) < 1) std::cout << "Inaccuracy: ";
            else if (reason == R_MISTAKE && fabs(delta) < 4) std::cout << "Mistake: ";
            else if (reason == R_MISTAKE && fabs(delta) < 10) std::cout << "Blunder: ";
            else if (reason == R_MISTAKE && fabs(delta) >= 10) std::cout << "Massive blunder: ";
            else if (reason == R_NONE) std::cout << "(Warning) No reason: ";
            else std::cout << "(Error) Invalid reason " << reason << ": ";
            std::cout << delta << std::endl;
            if (reason == R_MISTAKE && bestMvName != "") std::cout << "Best move was: " << bestMvName << std::endl;
        }
        scoreByReason[reason] += delta;
        expScore = newExpScore;
    }

    void printScoreByReason()
    {
        std::cout << "Baseline score: " << getInitialScore() << std::endl;
        std::cout << "Score due to luck: " << scoreByReason[R_LUCK] << std::endl;
        std::cout << "Score due to mistakes: " << scoreByReason[R_MISTAKE] << std::endl;
        if (fabs(scoreByReason[R_NONE]) > EPS) std::cout << "(Warning) Score for no reason: " << scoreByReason[R_NONE] << std::endl;
    }
};

void printOccurs(const Occurs& occurs)
{
    bool first = true;
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        for (int j = 0; j < occurs[i]; ++j)
        {
            if (!first) std::cout << " ";
            std::cout<< i + 1;
            first = false;
        }
    }
    std::cout << std::endl;
}

Move chooseMove()
{
    std::cout << "Move: ";
    std::string s = "";
    while (s == "")
    {
        std::getline(std::cin, s);
    }
    std::stringstream ss(stringToLower(s));

    std::string name = "";
    std::vector<int> args;
    while (ss)
    {
        std::string word;
        ss >> word;
        if (word == "") continue;
        if (isNumber(word)) args.push_back(stoi(word));
        else if (name == "") name = word;
        else name += " " + word;
    }

    return Move(name, args);
}

const std::string numNames[NUM_SIDES] = {"ones", "twos", "threes", "fours", "fives", "sixes"};
void simCollMoves(State& s, int num, int collected, int verbosity, bool printEvalLM, bool manualRolls, bool manualMoves)
{
    bool first = true;
    bool cont = true;
    while (cont)
    {
        cont = false;
        int left = NUM_DICE - collected;
        Move bestMv = getMoveColl(s.free, s.goods, num, left);

        s.updateExpScore(num * collected + bestMv.score, first ? R_NONE : R_LUCK, printEvalLM);
        first = false;

        if (verbosity >= 4)
        {
            std::cout << "Number of " << numNames[num - 1] << " collected: " << collected << std::endl;
            std::cout << "Goods: " << s.goods << std::endl;
        }

        moveSelectColl:

        Move mv = manualMoves ? chooseMove() : bestMv;

        if (!manualMoves && verbosity >= 3)
        {
            std::cout << "Move: " << mv.toString() << std::endl;
        }

        switch (mv.id)
        {
        case M_STOP_COLL:
            s.updateExpScore(num * collected + getScore(s.free, s.goods), R_MISTAKE, printEvalLM);
            break;
        
        case M_CONT_COLL:
            if (s.goods && left)
            {
                s.updateExpScore(num * collected + getScoreCollContinue(s.free, s.goods, num, left), R_MISTAKE, printEvalLM);
                cont = true;
                --s.goods;
                int newHits = 0;
                if (!manualRolls)
                {
                    for (int i = 0; i < left; ++i)
                    {
                        if (randDiceRoll() == 0) ++newHits;
                    }
                    if (verbosity >= 3)
                    {
                        std::cout << "Number of " << numNames[num - 1] << " rolled: " << newHits << std::endl;
                    }
                }
                else
                {
                    newHits = -1;
                    std::cout << "Number of " << numNames[num - 1] << " rolled: ";
                    while (newHits < 0 || newHits > left)
                    {
                        std::cin >> newHits;
                    }
                }
                collected += newHits;
            }
            else s.fail = true;
            break;

        default:
            s.fail = true;
        }

        if (manualMoves && s.fail)
        {
            std::cout << "Invalid move." << std::endl;
            s.fail = false;
            goto moveSelectColl;
        }
    }
    
    if (verbosity >= 4)
    {
        std::cout << "Number of " << numNames[num - 1] << " collected: " << collected << std::endl;
    }

    int reward = num * collected;
    s.score += reward;

    if (verbosity >= 2)
    {
        std::cout << "Won " << reward << " points." << std::endl;
    }
}

bool isDone(const Occurs& occurs)
{
    return std::all_of(occurs.begin(), occurs.end(), [](int n) {return n == 0;});
}

void simRegMove(State& s, Occurs& occurs, int ed, int reward, int verbosity, bool printEvalLM, bool manualRolls)
{
    bool instaDone = isDone(occurs);
    int rolls = 0;
    bool won;

    if (!manualRolls)
    {
        while (!isDone(occurs) && rolls < s.goods)
        {
            if (verbosity >= 4)
            {
                std::cout << "Need: ";
                printOccurs(occurs);
            }

            ++rolls;
            randRemOcc(ed, occurs);
        }
        won = isDone(occurs);
    }
    else
    {
        if (!isDone(occurs))
        {
            std::cout << "Number of rolls to complete the combination (-1 for if you didn't): ";
            std::cin >> rolls;
        }

        won = rolls >= 0 && rolls <= s.goods;
        if (!won) rolls = s.goods;
    }

    s.goods -= rolls;
    if (won) s.score += reward;

    if (!manualRolls && verbosity >= 2)
    {
        if (won) std::cout << "Took " << rolls << " rolls to complete the combination. " << std::endl;
        else std::cout << "Didn't complete the combination." << std::endl;
    }

    s.updateExpScore(getScore(s.free, s.goods), instaDone ? R_NONE : R_LUCK, printEvalLM);

    if (won && verbosity >= 3)
    {
        std::cout << "Won " << reward << " points." << std::endl;
    }
}

void chooseOcc(Occurs& occurs)
{
    occurs.fill(0);
    std::cout << "Rolled: ";
    int side;
    for (int i = 0; i < NUM_DICE; ++i)
    {
        std::cin >> side;
        if (side >= 1 && side <= NUM_SIDES) ++occurs[side - 1];
        else --i;
    }
}

int simGame(int verbosity=-1, bool printEvalLM=false, bool manualRolls=false, bool manualMoves=false)
{
    State s;
    Occurs diceOccurs;

    while (!s.fail && s.getTurn() <= NUM_COMBOS)
    {
        if (verbosity >= 3)
        {
            std::cout << "Turn: " << s.getTurn() << "/" << NUM_COMBOS << std::endl;
        }

        if (verbosity >= 1)
        {
            std::cout << "Current score: " << s.score << std::endl;
        }
        
        s.updateExpScore(getScore(s.free, s.goods), R_NONE, printEvalLM);

        if (printEvalLM)
        {
            std::cout << "Expected final score: " << s.expScore << std::endl;
        }

        s.goods += GOODS_PER_TURN;

        diceRoll:

        if (!manualRolls) randOcc(diceOccurs);
        else chooseOcc(diceOccurs);

        if (!manualRolls && verbosity >= 2)
        {
            std::cout << "Rolled: ";
            printOccurs(diceOccurs);
        }

        Move bestMv = getMove(s.free, s.goods, diceOccurs);

        s.updateExpScore(bestMv.score, R_LUCK, printEvalLM);

        if (verbosity >= 1)
        {
            std::cout << "Goods: " << s.goods << std::endl;
        }

        moveSelect:

        Move mv = manualMoves ? chooseMove() : bestMv;

        if (!manualMoves && verbosity >= 1)
        {
            std::cout << "Move: " << mv.toString() << std::endl;
        }

        int id = mv.id;
        switch (id)
        {
        case M_ERROR:
        case M_STOP_COLL:
        case M_CONT_COLL:
            s.fail = true;
            break;
        
        case M_REROLL:
            if (s.goods)
            {
                --s.goods;
                s.updateExpScore(getScoreCont(s.free, s.goods), R_MISTAKE, printEvalLM, bestMv.toString());
                goto diceRoll;
            }
            else s.fail = true;
            break;
        
        default:
            if (id >= 0 && id < NUM_COMBOS && isFree(s.free, id))
            {
                s.free = setUsed(s.free, id);
                int num = combos[id]->getCollectNumber();
                if (num)
                {
                    int collected = diceOccurs[num - 1];
                    s.updateExpScore(num * collected + getScoreColl(s.free, s.goods, num, NUM_DICE - collected), R_MISTAKE, printEvalLM, bestMv.toString());
                    simCollMoves(s, num, collected, verbosity, printEvalLM, manualRolls, manualMoves);
                }
                else
                {
                    Target t = combos[id]->getTargetByArgs(mv.args);
                    Occurs newOccurs = t.occurs;
                    int extraDice = calcExtraDice(newOccurs);
                    remOcc(diceOccurs, newOccurs);
                    s.updateExpScore(getMoveScore(s.free, s.goods, newOccurs, extraDice, t.points), R_MISTAKE, printEvalLM, bestMv.toString());
                    simRegMove(s, newOccurs, extraDice, t.points, verbosity, printEvalLM, manualRolls);
                }
            }
            else s.fail = true;
        }

        if (manualMoves && s.fail)
        {
            std::cout << "Invalid move" << std::endl;
            s.fail = false;
            goto moveSelect;
        }

        if (verbosity >= 1) std::cout << std::endl;
    }

    if (s.fail)
    {
        std::cout << "Encountered error!" << std::endl;
        return 0;
    }

    s.score += s.goods * POINTS_PER_GOOD;

    s.updateExpScore(0, R_NONE, printEvalLM);

    if (verbosity >= 0)
    {
        std::cout << "Final score: " << s.score << std::endl;
    }

    if (printEvalLM)
    {
        s.printScoreByReason();
    }

    return s.score;
}

struct Stats
{
    double mean;
    double stdev;
    int median;
    int perc5;
    int perc95;
    std::vector<int> modes;
};

Stats findStats(int n)
{
    std::vector<int> res;
    for (int i = 0; i < n; ++i)
    {
        res.push_back(simGame());
    }
    std::sort(res.begin(), res.end());

    Stats stats;

    long long sum = 0;
    long long sqSum = 0;
    int maxCnt = 0;
    int cnt = -1;
    int curr = -1;

    for (int r : res)
    {
        sum += r;
        sqSum += r * r;

        if (curr != r)
        {
            cnt = 0;
            curr = r;
        }

        ++cnt;
        if (cnt > maxCnt)
        {
            stats.modes.clear();
            maxCnt = cnt;
        }

        if (cnt == maxCnt) stats.modes.push_back(curr);
    }

    stats.mean = (double) sum / n;
    stats.stdev = sqrt((double) sqSum / n - stats.mean * stats.mean);
    stats.median = res[n / 2];
    stats.perc5 = res[n / 20];
    stats.perc95 = res[19 * n / 20];

    return stats;
}

void storeModel()
{
    std::string name = "model_" + std::to_string(NUM_TRIALS_2);
    std::ofstream file(name.c_str());

    std::cout << "Storing the model in " << name << "." << std::endl;

    for (int i = 0; i < DIFF_CODES; ++i)
    {
        file << indexToCode[i] << ' ';
    }

    for (int i = 0; i <= MAX_EXTRA_DICE; ++i)
    {
        for (int j = 0; j < DIFF_CODES; ++j)
        {
            for (int k = 0; k <= MAX_GOODS; ++k)
            {
                file << rollsDistr[i][j][k] << ' ';
            }
        }
    }

    for (int i = 0; i <= NUM_DICE; ++i)
    {
        for (int j = 0; j <= NUM_DICE; ++j)
        {
            file << leftDistr[i][j] << ' ';
        }
    }

    for (int i = 0; i < NUM_MASKS; ++i)
    {
        for (int j = 0; j <= MAX_GOODS; ++j)
        {
            file << isFound[i][j] << ' ';
            if (isFound[i][j]) file << score[i][j] << ' ';
        }
    }

    std::cout << "Stored the model." << std::endl;
}

bool loadModel()
{
    std::string name = "model_" + std::to_string(NUM_TRIALS_2);
    std::ifstream file(name.c_str());

    if (!file) return false;

    std::cout << "Loading the model from " << name << "." << std::endl;

    for (int i = 0; i < DIFF_CODES; ++i)
    {
        file >> indexToCode[i];
        codeToIndex[indexToCode[i]] = i;
    }

    for (int i = 0; i <= MAX_EXTRA_DICE; ++i)
    {
        for (int j = 0; j < DIFF_CODES; ++j)
        {
            for (int k = 0; k <= MAX_GOODS; ++k)
            {
                file >> rollsDistr[i][j][k];
            }
        }
    }

    for (int i = 0; i <= NUM_DICE; ++i)
    {
        for (int j = 0; j <= NUM_DICE; ++j)
        {
            file >> leftDistr[i][j];
        }
    }

    for (int i = 0; i < NUM_MASKS; ++i)
    {
        for (int j = 0; j <= MAX_GOODS; ++j)
        {
            file >> isFound[i][j];
            if (isFound[i][j]) file >> score[i][j];
        }
    }

    std::cout << "Loaded the model." << std::endl;
    return true;
}

void fitModel()
{
    std::cout << "Fitting the model." << std::endl;
    generateLefts();
    findRollsDistr();
    findLeftDistr();
    getInitialScore();
    std::cout << "\nConsidered " << statesVisCnt << " states." << std::endl;
}

const std::string help = "Possible commands: play / p, example / e, test, expected, credits, help, exit.";
const std::string credits = "Made by Emil Indzhev.";

void shell()
{
    std::string cmd;
    std::cout << "\n" << help << std::endl;

    while (true)
    {
        std::cout << "\nEnter command: ";
        std::cin >> cmd;

        if (cmd == "play" || cmd == "p")
        {
            bool manualRolls;
            bool manualMoves;
            bool printEvalLM;

            std::cout << "Manual rolls (0 / 1): ";
            std::cin >> manualRolls;

            std::cout << "Manual moves (0 / 1): ";
            std::cin >> manualMoves;

            std::cout << "Evaluate luck and mistakes (0 / 1): ";
            std::cin >> printEvalLM;

            std::cout << std::endl;
            simGame(3, printEvalLM, manualRolls, manualMoves);
        }
        else if (cmd == "example" || cmd == "e")
        {
            int verbosity;
            std::cout << "Verbosity (0 - 4): ";
            std::cin >> verbosity;
            std::cout << std::endl;
            simGame(verbosity);
        }
        else if (cmd == "test")
        {
            int numTests;
            std::cout << "Number of tests: ";
            std::cin >> numTests;
            Stats stats = findStats(numTests);
            std::cout << std::endl;
            std::cout << "Mean: " << stats.mean << std::endl;
            std::cout << "Stdev: " << stats.stdev << std::endl;
            std::cout << "5th  %ile: " << stats.perc5 << std::endl;
            std::cout << "Median: " << stats.median << std::endl;
            std::cout << "95th  %ile: " << stats.perc95 << std::endl;
            std::cout << "Modes:";
            for (int m : stats.modes)
            {
                std::cout << " " << m;
            }
            std::cout << std::endl;
        }
        else if (cmd == "expected")
        {
            std::cout << "Expected score: " << getScore(NUM_MASKS - 1, 0) << std::endl;
        }
        if (cmd == "credits")
        {
            std::cout << credits << std::endl;
        }
        if (cmd == "help")
        {
            std::cout << help << std::endl;
        }
        if (cmd == "exit") break;
    }
}

int main()
{
    generator.seed(time(nullptr));

    std::cout << "Number of trails per state, higher leads to a better model: ";
    std::cin >> NUM_TRIALS_2;

    if (!loadModel())
    {
        fitModel();
        storeModel();
    }

    shell();

    return 0;
}
