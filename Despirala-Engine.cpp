#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <random>
#include <chrono>
#include <math.h>
#include <time.h>

const double IO_EPS = 5e-4;
const int IO_PRECISION = 3;
const int FILE_PRECISION = 9;

const double INF = 1e6;

// Game configuration
const int NUM_DICE = 6;
const int NUM_SIDES = 6;
const int NUM_COMBOS = 14;
const int GOODS_PER_TURN = 5;
const int POINTS_PER_GOOD = 1;

// Depends on game configuration
const int MAX_CODE = 1 << (NUM_SIDES + NUM_DICE);
const int VALID_CODES = 30; // Sum_{0 <= k <= NUM_DICE} #ways to partition k into at most NUM_SIDES partitions
const int VALID_ORD_CODES = 924; // Sum_{0 <= k <= NUM_DICE} #ways to partition k into NUM_SIDES ordered partitions
const int STARTS_ORD_CODES[NUM_DICE + 2] = {0, 1, 7, 28, 84, 210, 462, 924}; // 0 <= t <= NUM_DICE + 1: Sum_{0 <= k <= t} #ways to partition k into NUM_SIDES ordered partitions
const int NUM_MASKS = 1 << NUM_COMBOS;
const int MAX_GOODS = NUM_COMBOS * GOODS_PER_TURN;

int NUM_TRIALS; // Trials per state

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
Combo* const combos[] = {new CollectCombo("Collect", 1),
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

int occursToCode(Occurs occurs, bool ordered)
{
    if (!ordered) std::sort(occurs.begin(), occurs.end(), std::greater<int>());

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

bool codeToOccurs(int code, Occurs& occurs, bool ordered)
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
        if (!ordered && occurs[i - 1] < occurs[i]) return false;
    }
    return true;
}

std::mt19937 generator;
std::uniform_int_distribution<int> diceDistribution(0, NUM_SIDES - 1);

int randDiceRoll()
{
    return diceDistribution(generator);
}

void randOcc(Occurs& occurs, int numDice = NUM_DICE)
{
    occurs.fill(0);
    for (int i = 0; i < numDice; ++i)
    {
        ++occurs[randDiceRoll()];
    }
}

bool remOcc(Occurs& occurs, const Occurs& diceOccurs)
{
    bool changed = false;
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        if (occurs[i] > 0 && diceOccurs[i] > 0) changed = true;
        occurs[i] = std::max(occurs[i] - diceOccurs[i], 0);
    }
    return changed;
}

int diceInOcc(const Occurs& occurs)
{
    return std::accumulate(occurs.begin(), occurs.end(), 0);
}

bool randRemOcc(int extraDice, Occurs& occurs)
{
    Occurs diceOccurs;
    int numDice = extraDice + diceInOcc(occurs);
    randOcc(diceOccurs, numDice);
    return remOcc(occurs, diceOccurs);
}

std::unordered_map<int, int> codeToIndex;
int indexToCode[VALID_CODES];
int diceInIndex[VALID_CODES];

std::unordered_map<int, int> ordCodeToIndex;
int ordIndexToCode[VALID_ORD_CODES];
double ordIndexProb[VALID_ORD_CODES];

int fact(int n)
{
    return n ? n * fact(n - 1) : 1;
}

double occursProb(const Occurs& occurs)
{
    int numDice = diceInOcc(occurs);
    double prob = fact(numDice);
    for (int cnt : occurs)
    {
        prob /= fact(cnt);
    }
    return prob / pow(NUM_SIDES, numDice);
}

void genCodeIdxMap()
{
    Occurs occurs;
    int lastIdx = 0;
    int ordLastIdx = 0;
    for (int i = 0; i <= MAX_CODE; ++i)
    {
        if (codeToOccurs(i, occurs, false))
        {
            codeToIndex[i] = lastIdx;
            indexToCode[lastIdx] = i;
            diceInIndex[lastIdx++] = diceInOcc(occurs);
        }

        if (codeToOccurs(i, occurs, true))
        {
            ordCodeToIndex[i] = ordLastIdx;
            ordIndexToCode[ordLastIdx] = i;
            ordIndexProb[ordLastIdx++] = occursProb(occurs);
        }
    }
}

double rollsDistr[NUM_DICE + 1][VALID_CODES][MAX_GOODS + 1];

void findRollsDistrSingle(int extraDice, int idx)
{
    Occurs occurs;
    double changeProb = 0;
    int code = indexToCode[idx];
    int numDice = extraDice + diceInIndex[idx];
    codeToOccurs(code, occurs, false);
    for (int i = STARTS_ORD_CODES[numDice]; i < STARTS_ORD_CODES[numDice + 1]; ++i)
    {
        Occurs diceOccurs;
        int diceCode = ordIndexToCode[i];
        codeToOccurs(diceCode, diceOccurs, true);
        Occurs newOccurs = occurs;
        if (!remOcc(newOccurs, diceOccurs)) continue;
        changeProb += ordIndexProb[i];
        int newCode = occursToCode(newOccurs, false);
        int newIdx = codeToIndex[newCode];
        for (int j = 0; j <= MAX_GOODS; ++j)
        {
            rollsDistr[extraDice][idx][j] += rollsDistr[extraDice][newIdx][j] * ordIndexProb[i];
        }
    }
    for (int i = MAX_GOODS; i >= 0 ; --i)
    {
        rollsDistr[extraDice][idx][i] = 0;
        for (int j = i - 1; j >= 0 ; --j)
        {
            rollsDistr[extraDice][idx][i] += rollsDistr[extraDice][idx][j] * pow(1 - changeProb, i - j - 1);
        }
    }
}

void findRollsDistr()
{
    for (int extraDice = 0; extraDice <= NUM_DICE; ++extraDice)
    {
        rollsDistr[extraDice][0][0] = 1;
        for (int i = 1; i < VALID_CODES; ++i)
        {
            if (extraDice + diceInIndex[i] > NUM_DICE) continue; 
            findRollsDistrSingle(extraDice, i);
        }
    }
}

double leftDistr[NUM_DICE + 1][NUM_DICE + 1];

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

double getMoveScore(int free, int goods, const Occurs& occurs, int extraDice, int reward)
{
    double score = 0;
    int code = occursToCode(occurs, false);
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
                int extraDice = NUM_DICE - diceInOcc(newOccurs);
                remOcc(newOccurs, diceOccurs);
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
    else if (NUM_TRIALS > 0)
    {
        for (int t = 0; t < NUM_TRIALS; ++t)
        {
            Occurs diceOccurs;
            randOcc(diceOccurs);
            Move mov = getMove(free, goods, diceOccurs);
            score[free][goods] += mov.score / NUM_TRIALS;
        }
    }
    else
    {
        for (int i = STARTS_ORD_CODES[NUM_DICE]; i < STARTS_ORD_CODES[NUM_DICE + 1]; ++i)
        {
            Occurs diceOccurs;
            int diceCode = ordIndexToCode[i];
            codeToOccurs(diceCode, diceOccurs, true);
            Move mov = getMove(free, goods, diceOccurs);
            score[free][goods] += mov.score * ordIndexProb[i];
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

    void updateExpScore(double newExpScore, int reason, bool evalLM, const std::string& bestMvName = "")
    {
        newExpScore += score;
        double delta = newExpScore - expScore;
        bool print = evalLM && (reason == R_LUCK || fabs(newExpScore - expScore) > IO_EPS);
        if (print)
        {
            if (reason == R_LUCK && delta > IO_EPS) std::cout << "Good luck: ";
            else if (reason == R_LUCK && delta < -IO_EPS) std::cout << "Bad luck: ";
            else if (reason == R_LUCK && fabs(delta) <= IO_EPS) std::cout << "Neutral luck: ";
            else if (reason == R_MISTAKE && -delta < 1) std::cout << "Inaccuracy: ";
            else if (reason == R_MISTAKE && -delta < 4) std::cout << "Mistake: ";
            else if (reason == R_MISTAKE && -delta < 10) std::cout << "Blunder: ";
            else if (reason == R_MISTAKE && -delta >= 10) std::cout << "Massive blunder: ";
            else if (reason == R_NONE) std::cout << "(Warning) No reason: ";
            else std::cout << "(Error) Invalid reason " << reason << ": ";
            std::cout << std::fixed << std::setprecision(IO_PRECISION) << delta << std::endl;
            if (reason == R_MISTAKE && bestMvName != "") std::cout << "Best move was: " << bestMvName << std::endl;
        }
        scoreByReason[reason] += delta;
        expScore = newExpScore;
    }

    void printScoreByReason()
    {
        std::cout << "Baseline score: " << std::fixed << std::setprecision(IO_PRECISION) << getInitialScore() << std::endl;
        std::cout << "Score due to luck: " << std::fixed << std::setprecision(IO_PRECISION) << scoreByReason[R_LUCK] << std::endl;
        std::cout << "Score due to mistakes: " << std::fixed << std::setprecision(IO_PRECISION) << scoreByReason[R_MISTAKE] << std::endl;
        if (fabs(scoreByReason[R_NONE]) > IO_EPS) std::cout << "(Warning) Score for no reason: " << std::fixed << std::setprecision(IO_PRECISION) << scoreByReason[R_NONE] << std::endl;
    }
};

#define LOG_NONE 0
#define LOG_READ 1
#define LOG_WRITE 2

struct Config
{
    bool verbose;
    bool manualRolls;
    bool manualMoves;
    bool evalLM;
    int logMode;
    std::ifstream logIn;
    std::ofstream logOut;

    Config(bool verbose = false, int logMode = LOG_NONE, const std::string logName = "", bool evalLM = false, bool manualRolls = false, bool manualMoves = false):
        verbose(verbose),
        manualRolls(manualRolls),
        manualMoves(manualMoves),
        evalLM(evalLM),
        logMode(logMode)
    {
        if (logMode == LOG_READ) logIn.open(logName.c_str());
        else if (logMode == LOG_WRITE) logOut.open(logName.c_str());
    }  
};

void printOcc(const Occurs& occurs, Config& config, bool useConfig)
{
    if (!useConfig || config.logMode != LOG_WRITE) std::cout << "Rolled: ";
    std::ostream& out = useConfig && config.logMode == LOG_WRITE ? config.logOut : std::cout;

    bool first = true;
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        for (int j = 0; j < occurs[i]; ++j)
        {
            if (!first) out << " ";
            out<< i + 1;
            first = false;
        }
    }
    out << std::endl;
}

void chooseOcc(Occurs& occurs, Config& config)
{
    if (config.logMode != LOG_READ) std::cout << "Rolled: ";
    std::istream& in = config.logMode == LOG_READ ? config.logIn : std::cin;

    occurs.fill(0);
    int side;
    for (int i = 0; i < NUM_DICE; ++i)
    {
        in >> side;
        if (side >= 1 && side <= NUM_SIDES) ++occurs[side - 1];
        else --i;
    }
}

Move chooseMove(Config& config)
{
    if (config.logMode != LOG_READ) std::cout << "Move: ";
    std::istream& in = config.logMode == LOG_READ ? config.logIn : std::cin;

    std::string s = "";
    while (s == "")
    {
        std::getline(in, s);
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

int chooseNumRolls(Config& config)
{
    if (config.logMode != LOG_READ) std::cout << "Number of rolls to complete the combination (-1 for if you didn't): ";
    std::istream& in = config.logMode == LOG_READ ? config.logIn : std::cin;

    int rolls;
    in >> rolls;
    return rolls;
}

const std::string numNames[NUM_SIDES] = {"ones", "twos", "threes", "fours", "fives", "sixes"};
int chooseNumNewHits(int num, int left, Config& config)
{
    if (config.logMode != LOG_READ) std::cout << "Number of " << numNames[num - 1] << " rolled: ";
    std::istream& in = config.logMode == LOG_READ ? config.logIn : std::cin;

    int newHits = -1;
    while (newHits < 0 || newHits > left)
    {
        in >> newHits;
    }
    return newHits;
}

void simCollMoves(State& s, int num, int collected, Config& config)
{
    bool first = true;
    bool cont = true;
    while (cont)
    {
        cont = false;
        int left = NUM_DICE - collected;
        Move bestMv = getMoveColl(s.free, s.goods, num, left);

        s.updateExpScore(num * collected + bestMv.score, first ? R_NONE : R_LUCK, config.evalLM);
        first = false;

        moveSelectColl:

        Move mv = config.manualMoves || config.logMode == LOG_READ ? chooseMove(config) : bestMv;

        if (!config.manualMoves && config.verbose) std::cout << "Move: " << mv.toString() << std::endl;

        switch (mv.id)
        {
        case M_STOP_COLL:
            s.updateExpScore(num * collected + getScore(s.free, s.goods), R_MISTAKE, config.evalLM);
            break;
        
        case M_CONT_COLL:
            if (s.goods && left)
            {
                s.updateExpScore(num * collected + getScoreCollContinue(s.free, s.goods, num, left), R_MISTAKE, config.evalLM);
                cont = true;
                --s.goods;
                int newHits = 0;
                if (config.manualRolls || config.logMode == LOG_READ) newHits = chooseNumNewHits(num, left, config);
                else
                {
                    for (int i = 0; i < left; ++i)
                    {
                        if (randDiceRoll() == 0) ++newHits;
                    }
                }

                if (config.logMode == LOG_WRITE) config.logOut << newHits << std::endl;

                if (!config.manualRolls && config.verbose)
                {
                    std::cout << "Number of " << numNames[num - 1] << " rolled: " << newHits << std::endl;
                }

                collected += newHits;
            }
            else s.fail = true;
            break;

        default:
            s.fail = true;
        }

        if (config.logMode == LOG_WRITE) config.logOut << mv.toString() << std::endl;

        if (config.manualMoves && s.fail)
        {
            std::cout << "Invalid move." << std::endl;
            s.fail = false;
            goto moveSelectColl;
        }
    }

    int reward = num * collected;
    s.score += reward;

    if (config.verbose)
    {
        std::cout << "Won " << reward << " points." << std::endl;
    }
}

bool isDone(const Occurs& occurs)
{
    return std::all_of(occurs.begin(), occurs.end(), [](int n) {return n == 0;});
}

void simRegMove(State& s, Occurs& occurs, int extraDice, int reward, Config& config)
{
    bool instaDone = isDone(occurs);
    int rolls = 0;
    bool won;

    if (config.manualRolls || config.logMode == LOG_READ)
    {
        if (!isDone(occurs)) rolls = chooseNumRolls(config);

        won = rolls >= 0 && rolls <= s.goods;
        if (!won) rolls = s.goods;
    }
    else
    {
        while (!isDone(occurs) && rolls < s.goods)
        {
            ++rolls;
            randRemOcc(extraDice, occurs);
        }
        won = isDone(occurs);
    }

    if (!instaDone && config.logMode == LOG_WRITE) config.logOut << (won ? rolls : -1) << std::endl;

    if (!config.manualRolls && config.verbose)
    {
        if (won) std::cout << "Took " << rolls << " rolls to complete the combination. " << std::endl;
        else std::cout << "Didn't complete the combination." << std::endl;
    }

    if (won && config.verbose)
    {
        std::cout << "Won " << reward << " points." << std::endl;
    }

    s.goods -= rolls;
    if (won) s.score += reward;

    s.updateExpScore(getScore(s.free, s.goods), instaDone ? R_NONE : R_LUCK, config.evalLM);
}

int simGame(Config& config)
{
    State s;
    Occurs diceOccurs;

    if (config.verbose) std::cout << std::endl;

    while (!s.fail && s.getTurn() <= NUM_COMBOS)
    {
        if (config.verbose)
        {
            std::cout << "Turn: " << s.getTurn() << "/" << NUM_COMBOS << std::endl;
            std::cout << "Current score: " << s.score << std::endl;
        }

        s.updateExpScore(getScore(s.free, s.goods), R_NONE, config.evalLM);

        if (config.evalLM)
        {
            std::cout << "Expected final score: " << std::fixed << std::setprecision(IO_PRECISION) << s.expScore << std::endl;
        }

        s.goods += GOODS_PER_TURN;

        diceRoll:

        if (config.manualRolls || config.logMode == LOG_READ) chooseOcc(diceOccurs, config);
        else randOcc(diceOccurs);

        if (config.logMode == LOG_WRITE) printOcc(diceOccurs, config, true);
        if (!config.manualRolls && config.verbose) printOcc(diceOccurs, config, false);

        Move bestMv = getMove(s.free, s.goods, diceOccurs);

        s.updateExpScore(bestMv.score, R_LUCK, config.evalLM);

        if (config.verbose)
        {
            std::cout << "Goods: " << s.goods << std::endl;
        }

        moveSelect:

        Move mv = config.manualMoves || config.logMode == LOG_READ ? chooseMove(config) : bestMv;

        if (!config.manualMoves && config.verbose) std::cout << "Move: " << mv.toString() << std::endl;

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
                if (config.logMode == LOG_WRITE) config.logOut << mv.toString() << std::endl;

                --s.goods;
                s.updateExpScore(getScoreCont(s.free, s.goods), R_MISTAKE, config.evalLM, bestMv.toString());
                goto diceRoll;
            }
            else s.fail = true;
            break;
        
        default:
            if (id >= 0 && id < NUM_COMBOS && isFree(s.free, id))
            {
                if (config.logMode == LOG_WRITE) config.logOut << mv.toString() << std::endl;

                s.free = setUsed(s.free, id);
                int num = combos[id]->getCollectNumber();
                if (num)
                {
                    int collected = diceOccurs[num - 1];
                    s.updateExpScore(num * collected + getScoreColl(s.free, s.goods, num, NUM_DICE - collected), R_MISTAKE, config.evalLM, bestMv.toString());
                    simCollMoves(s, num, collected, config);
                }
                else
                {
                    Target t = combos[id]->getTargetByArgs(mv.args);
                    Occurs newOccurs = t.occurs;
                    int extraDice = NUM_DICE - diceInOcc(newOccurs);
                    remOcc(newOccurs, diceOccurs);
                    s.updateExpScore(getMoveScore(s.free, s.goods, newOccurs, extraDice, t.points), R_MISTAKE, config.evalLM, bestMv.toString());
                    simRegMove(s, newOccurs, extraDice, t.points, config);
                }
            }
            else s.fail = true;
        }

        if (config.manualMoves && s.fail)
        {
            std::cout << "Invalid move" << std::endl;
            s.fail = false;
            goto moveSelect;
        }

        if (config.verbose) std::cout << std::endl;
    }

    if (s.fail)
    {
        std::cout << "Encountered error!" << std::endl;
        return 0;
    }

    s.score += s.goods * POINTS_PER_GOOD;

    s.updateExpScore(0, R_NONE, config.evalLM);

    if (config.verbose)
    {
        std::cout << "Final score: " << s.score << std::endl;
    }

    if (config.evalLM)
    {
        s.printScoreByReason();
    }

    return s.score;
}

struct Stats
{
    double mean;
    double stdev;
    int perc5;
    int perc25;
    int perc50;
    int perc75;
    int perc95;
    std::vector<int> modes;
    std::vector<int> distr;
};

Stats findStats(int n)
{
    Stats stats;

    for (int i = 0; i < n; ++i)
    {
        Config config;
        int r = simGame(config);

        while (r >= (int) stats.distr.size()) stats.distr.push_back(0);
        ++stats.distr[r];
    }

    long long sum = 0;
    long long sqSum = 0;
    int modeCnt = 0;
    int cumDistr = 0;

    for (int r = 0; r < (int) stats.distr.size(); ++r)
    {
        sum += (long long) stats.distr[r] * r;
        sqSum += (long long) stats.distr[r] * r * r;

        if (stats.distr[r] > modeCnt)
        {
            stats.modes.clear();
            modeCnt = stats.distr[r];
        }

        if (stats.distr[r] == modeCnt) stats.modes.push_back(r);

        if (cumDistr < n / 20) stats.perc5 = r;
        if (cumDistr < n / 4) stats.perc25 = r;
        if (cumDistr < n / 2) stats.perc50 = r;
        if (cumDistr < 3 * n / 4) stats.perc75 = r;
        if (cumDistr < 19 * n / 20) stats.perc95 = r;

        cumDistr += stats.distr[r];
    }

    stats.mean = (double) sum / n;
    stats.stdev = sqrt((double) sqSum / n - stats.mean * stats.mean);

    return stats;
}

std::string modelName()
{
    return (NUM_TRIALS > 0 ? std::to_string(NUM_TRIALS) : "max") + ".model";
}

void storeModel()
{
    std::string name = modelName();
    std::ofstream file(name.c_str());

    std::cout << "Storing model in " << name << "." << std::endl;

    for (int i = 0; i < VALID_CODES; ++i)
    {
        file << indexToCode[i] << ' ' << diceInIndex[i] << ' ';
    }

    for (int i = 0; i < VALID_ORD_CODES; ++i)
    {
        file << ordIndexToCode[i] << ' ' << std::setprecision(FILE_PRECISION) << ordIndexProb[i] << ' ';
    }

    for (int i = 0; i <= NUM_DICE; ++i)
    {
        for (int j = 0; j < VALID_CODES; ++j)
        {
            if (i + diceInIndex[j] > NUM_DICE) continue;
            for (int k = 0; k <= MAX_GOODS; ++k)
            {
                file << std::setprecision(FILE_PRECISION) << rollsDistr[i][j][k] << ' ';
            }
        }
    }

    for (int i = 0; i <= NUM_DICE; ++i)
    {
        for (int j = 0; j <= NUM_DICE; ++j)
        {
            file << std::setprecision(FILE_PRECISION) << leftDistr[i][j] << ' ';
        }
    }

    for (int i = 0; i < NUM_MASKS; ++i)
    {
        for (int j = 0; j <= MAX_GOODS; ++j)
        {
            file << isFound[i][j] << ' ';
            if (isFound[i][j]) file << std::setprecision(FILE_PRECISION) << score[i][j] << ' ';
        }
    }

    std::cout << "Stored model." << std::endl;
}

bool loadModel()
{
    std::string name = modelName();
    std::ifstream file(name.c_str());

    if (!file) return false;

    std::cout << "Loading model from " << name << "." << std::endl;

    for (int i = 0; i < VALID_CODES; ++i)
    {
        file >> indexToCode[i] >> diceInIndex[i];
        codeToIndex[indexToCode[i]] = i;
    }

    for (int i = 0; i < VALID_ORD_CODES; ++i)
    {
        file >> ordIndexToCode[i] >> ordIndexProb[i];
        ordCodeToIndex[ordIndexToCode[i]] = i;
    }

    for (int i = 0; i <= NUM_DICE; ++i)
    {
        for (int j = 0; j < VALID_CODES; ++j)
        {
            if (i + diceInIndex[j] > NUM_DICE) continue;
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

    std::cout << "Loaded model." << std::endl;
    return true;
}

void computeModel()
{
    std::cout << "Computing model." << std::endl;
    genCodeIdxMap();
    findRollsDistr();
    findLeftDistr();
    getInitialScore();
    std::cout << std::endl << "Considered " << statesVisCnt << " states." << std::endl;
}

long long timestamp()
{
    return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

const std::string help = "Possible commands: play, replay, test, expected, credits, help, exit.";
const std::string credits = "Made by Emil Indzhev.";

void shell()
{
    std::string cmd;
    std::string lastLogName;

    std::cout << std::endl << help << std::endl;

    while (true)
    {
        std::cout << std::endl << "Enter command: ";
        std::cin >> cmd;

        if (cmd == "play")
        {
            bool manualRolls;
            bool manualMoves;
            bool evalLM;

            std::cout << "Manual rolls (0 / 1): ";
            std::cin >> manualRolls;

            std::cout << "Manual moves (0 / 1): ";
            std::cin >> manualMoves;

            std::cout << "Evaluate luck and mistakes (0 / 1): ";
            std::cin >> evalLM;

            lastLogName = std::to_string(timestamp()) + ".log";

            Config config(true, LOG_WRITE, lastLogName, evalLM, manualRolls, manualMoves);
            simGame(config);
        }
        else if (cmd == "replay")
        {
            bool evalLM;
            std::string logName;

            std::cout << "Log name or 'last' for last one: ";
            std::cin >> logName;

            if (logName == "last") logName = lastLogName;

            std::cout << "Evaluate luck and mistakes (0 / 1): ";
            std::cin >> evalLM;

            Config config(true, LOG_READ, logName, evalLM);
            simGame(config);
        }
        else if (cmd == "test")
        {
            int numTests;
            bool exportToFile;
            std::string fileName;

            std::cout << "Number of tests: ";
            std::cin >> numTests;

            std::cout << "Export raw data to file (0 / 1): ";
            std::cin >> exportToFile;

            if (exportToFile)
            {    
                std::cout << "File name: ";
                std::cin >> fileName;
            }

            Stats stats = findStats(numTests);

            std::cout << std::endl;
            std::cout << "Mean: " << std::fixed << std::setprecision(IO_PRECISION) << stats.mean << std::endl;
            std::cout << "Stdev: " << std::fixed << std::setprecision(IO_PRECISION) << stats.stdev << std::endl;
            std::cout << "5th percentile: " << stats.perc5 << std::endl;
            std::cout << "25th percentile: " << stats.perc25 << std::endl;
            std::cout << "50th percentile: " << stats.perc50 << std::endl;
            std::cout << "75th percentile: " << stats.perc75 << std::endl;
            std::cout << "95th percentile: " << stats.perc95 << std::endl;
            std::cout << (stats.modes.size() == 1 ? "Mode:" : "Modes:");
            for (int m : stats.modes)
            {
                std::cout << " " << m;
            }
            std::cout << std::endl;

            if (exportToFile)
            {    
                std::ofstream file(fileName.c_str());
                for (int d : stats.distr)
                {
                    file << d << std::endl;
                }
            }
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

    std::cout << "Number of trails per state or 0 for exact model (which costs " << STARTS_ORD_CODES[NUM_DICE + 1] - STARTS_ORD_CODES[NUM_DICE] << " trials per state): ";
    std::cin >> NUM_TRIALS;

    if (!loadModel())
    {
        computeModel();
        storeModel();
    }

    shell();

    return 0;
}
