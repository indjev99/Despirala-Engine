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
#include <assert.h>

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
bool MISERE; // Normal play or misere play

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
        auto it = cache.find(code);
        if (it == cache.end())
        {
            it = cache.insert({code, rawGetTargets(diceOccurs)}).first;
        }
        return it->second;
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
        if (MISERE) std::reverse(diceOccursPairs.begin(), diceOccursPairs.end());

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
        std::sort(args.begin(), args.end());
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
        if (MISERE) std::reverse(diceOccursPairs.begin(), diceOccursPairs.end());

        Occurs currTemplate(templateOccurs);

        int bestPoints = !MISERE ? -INF : INF;
        std::vector<Target> ts;
        Occurs occurs;
        do
        {
            std::vector<int> args = makeOccurs(currTemplate, diceOccursPairs, occurs);
            int currPoints = evalOccurs(occurs);
            if (!MISERE ? currPoints > bestPoints : currPoints < bestPoints)
            {
                bestPoints = currPoints;
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

    Move(int id, const std::vector<int>& args):
        Move(id, args, !MISERE ? -INF : INF) {}

    Move(int id, int arg, double score):
        Move(id, std::vector<int>(1, arg), score) {}

    Move(int id, int arg):
        Move(id, std::vector<int>(1, arg)) {}

    Move(int id, double score):
        Move(id, std::vector<int>(0), score) {}

    Move(int id):
        Move(id, std::vector<int>(0)) {}

    Move():
        Move(M_ERROR) {}

    Move(const std::string& name, const std::vector<int>& args);
    
    std::string toString() const;
};

bool operator<(const Move& a, const Move& b)
{
    return !MISERE ? a.score < b.score : a.score > b.score;
}

bool isFree(int free, int idx)
{
    return (free >> idx) & 1;
}

int setUsed(int free, int idx)
{
    return free - (1 << idx);
}

std::vector<Move> getMoveOptions(int free, int goods, const Occurs& diceOccurs)
{
    std::vector<Move> options;
    for (int i = 0; i < NUM_COMBOS; ++i)
    {
        if (!isFree(free, i)) continue;
        int num = combos[i]->getCollectNumber();
        if (num > 0)
        {
            options.emplace_back(i, num);
        }
        else
        {
            std::vector<Target> ts = combos[i]->getTargets(diceOccurs);
            for (const Target& t : ts)
            {
                options.emplace_back(i, t.args);
            }
        }
    }
    if (goods > 0 && !MISERE)
    {
        options.emplace_back(M_REROLL);
    }
    return options;
}

std::vector<Move> getCollMoveOptions(int goods, int left)
{
    std::vector<Move> options;
    options.emplace_back(M_STOP_COLL);
    if (goods > 0 && left > 0) options.emplace_back(M_CONT_COLL);
    return options;
}

bool isFound[NUM_MASKS][MAX_GOODS + 1];
double score[NUM_MASKS][MAX_GOODS + 1];
bool isFoundColl[NUM_MASKS][MAX_GOODS + 1][NUM_SIDES][NUM_DICE + 1];
double collScore[NUM_MASKS][MAX_GOODS + 1][NUM_SIDES][NUM_DICE + 1];

double getScore(int free, int goods);
double getContScore(int free, int goods);
double getCollScore(int free, int goods, int num, int left);

double getCollContScore(int free, int goods, int num, int left)
{
    double sc = 0;
    for (int i = 0 ; i <= left; ++i)
    {
        sc += leftDistr[left][i] * ((left - i) * num + getCollScore(free, goods - 1, num, i));
    }
    return sc;
}

std::vector<Move> getCollMoveOptionsScored(int free, int goods, int num, int left)
{
    std::vector<Move> options = getCollMoveOptions(goods, left);
    for (Move& option : options)
    {
        if (option.id == M_STOP_COLL)
        {
            option.score = getScore(free, goods);
        }
        else if (option.id == M_CONT_COLL)
        {
            option.score = getCollContScore(free, goods, num, left);
        }
    }
    std::sort(options.rbegin(), options.rend());
    return options;
}

// Optimized version of getCollMoveOptionsScored(...)[0]
Move getCollMove(int free, int goods, int num, int left)
{
    Move best(M_STOP_COLL, getScore(free, goods));
    if (goods > 0 && left > 0)
    {
        Move option(M_CONT_COLL, getCollContScore(free, goods, num, left));
        best = std::max(best, option);
    }
    return best;
}

double getCollScore(int free, int goods, int num, int left)
{
    if (!isFoundColl[free][goods][num][left])
    {
        isFoundColl[free][goods][num][left] = true;
        collScore[free][goods][num][left] = getCollMove(free, goods, num, left).score;
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

std::vector<Move> getMoveOptionsScored(int free, int goods, const Occurs& diceOccurs)
{
    std::vector<Move> options = getMoveOptions(free, goods, diceOccurs);
    for (Move& option : options)
    {
        int id = option.id;
        int newFree = id != M_REROLL ? setUsed(free, id) : free;
        int num = id != M_REROLL ? combos[id]->getCollectNumber() : 0;

        if (id == M_REROLL)
        {
            option.score = getContScore(newFree, goods - 1);
        }
        else if (num > 0)
        {
            int curr = diceOccurs[num - 1];
            option.score = num * curr + getCollScore(newFree, goods, num, NUM_DICE - curr);
        }
        else
        {
            Target t = combos[id]->getTargetByArgs(option.args);
            Occurs newOccurs = t.occurs;
            int extraDice = NUM_DICE - diceInOcc(newOccurs);
            remOcc(newOccurs, diceOccurs);
            option.score = getMoveScore(newFree, goods, newOccurs, extraDice, t.points);
        }
    }
    std::sort(options.rbegin(), options.rend());
    return options;
}

// Optimized version of getMoveOptionsScored(...)[0]
Move getMove(int free, int goods, const Occurs& diceOccurs)
{
    Move best;
    for (int i = 0; i < NUM_COMBOS; ++i)
    {
        if (!isFree(free, i)) continue;
        int newFree = setUsed(free, i);
        int num = combos[i]->getCollectNumber();
        if (num > 0)
        {
            int curr = diceOccurs[num - 1];
            Move option(i, num, num * curr + getCollScore(newFree, goods, num, NUM_DICE - curr));
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
                Move option(i, t.args, getMoveScore(newFree, goods, newOccurs, extraDice, t.points));
                best = std::max(best, option);
            }
        }
    }
    if (goods > 0 && !MISERE)
    {
        Move option(M_REROLL, getContScore(free, goods - 1));
        best = std::max(best, option);
    }
    return best;
}

int statesVisCnt = 0;

double getContScore(int free, int goods)
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
    if (statesVisCnt % 10000 == 0) std::cout << "." << std::flush;
    return score[free][goods];
}

double getScore(int free, int goods)
{
    return free ? getContScore(free, goods + GOODS_PER_TURN) : getContScore(free, goods);
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
    Move(M_ERROR, args)
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

    double scoreByReason[NUM_REASONS];

    std::vector<int> othersCumDistr;

    State():
        free(NUM_MASKS - 1),
        goods(0),
        score(0),
        expScore(getInitialScore())
    {
        std::fill(scoreByReason, scoreByReason + NUM_REASONS, 0);
    }

    int getTurn()
    {
        int turn = 0;
        for (int i = 0; i < NUM_COMBOS; ++i)
        {
            if ((free >> i & 1) == 0) ++turn;
        }
        return turn;
    }

    void updateExpScore(double newExpScore, int reason, const std::string& bestMvName = "")
    {
        newExpScore += score;
        int mult = !MISERE ? 1 : -1;
        double delta = newExpScore - expScore;
        bool print = reason == R_LUCK || fabs(delta) > IO_EPS;
        if (print)
        {
            if (reason == R_LUCK && mult * delta > IO_EPS) std::cout << "Good luck: ";
            else if (reason == R_LUCK && mult * delta < -IO_EPS) std::cout << "Bad luck: ";
            else if (reason == R_LUCK && fabs(delta) <= IO_EPS) std::cout << "Neutral luck: ";
            else if (reason == R_MISTAKE && -mult * delta < 1) std::cout << "Inaccuracy: ";
            else if (reason == R_MISTAKE && -mult * delta < 4) std::cout << "Mistake: ";
            else if (reason == R_MISTAKE && -mult * delta < 10) std::cout << "Blunder: ";
            else if (reason == R_MISTAKE && -mult * delta >= 10) std::cout << "Massive blunder: ";
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
#define LOG_WRITE 1
#define LOG_READ 2

struct Config
{
    int numPlayers;

    bool verbose;
    bool evalLM;

    int logMode;
    std::ifstream logIn;
    std::ofstream logOut;

    std::vector<bool> manualRollsAll;
    std::vector<bool> manualMovesAll;
    std::vector<int> competitiveAll;

    int player;
    bool manualRolls;
    bool manualMoves;
    int competitive;

    Config(int numPlayers = 1, bool verbose = false, bool evalLM = false, int logMode = LOG_NONE, const std::string logName = "",
           const std::vector<bool>& manualRollsAll = {}, const std::vector<bool>& manualMovesAll = {}, const std::vector<int>& competitiveAll = {}):
        numPlayers(numPlayers),
        verbose(verbose),
        evalLM(evalLM),
        logMode(logMode)
    {
        if (logMode == LOG_WRITE) logOut.open(logName.c_str());
        else if (logMode == LOG_READ) logIn.open(logName.c_str());

        if (logMode == LOG_READ) assert(numPlayers == 0);
        else assert(numPlayers > 0);

        if (logMode == LOG_WRITE) logOut << (this->numPlayers) << std::endl;
        else if (logMode == LOG_READ) logIn >> (this->numPlayers);

        this->manualRollsAll = vectorOrDefault(manualRollsAll);
        this->manualMovesAll = vectorOrDefault(manualMovesAll);
        this->competitiveAll = vectorOrDefault(competitiveAll);

        for (int player = 0; player < this->numPlayers; ++player)
        {
            assert(logMode != LOG_READ || !this->manualRollsAll[player]);
            assert(logMode != LOG_READ || !this->manualMovesAll[player]);
            assert(logMode != LOG_READ || !this->competitiveAll[player]);
            assert(!this->manualMovesAll[player] || !this->competitiveAll[player]);
        }

        setPlayer(0);
    }

    bool isDefault()
    {
        return numPlayers == 1 && !verbose && !evalLM && logMode == LOG_NONE && !manualRolls && !manualMoves;
    }

    void setPlayer(int player)
    {
        this->player = player;
        manualRolls = manualRollsAll[player];
        manualMoves = manualMovesAll[player];
        competitive = competitiveAll[player];
    }

private:

    std::vector<bool> vectorOrDefault(std::vector<bool> vec)
    {
        assert(vec.empty() || (int) vec.size() == numPlayers);

        if (vec.empty()) return std::vector<bool>(numPlayers, false);
        else return vec;
    }

    std::vector<int> vectorOrDefault(std::vector<int> vec)
    {
        assert(vec.empty() || (int) vec.size() == numPlayers);

        if (vec.empty()) return std::vector<int>(numPlayers, 0);
        else return vec;
    }
};

double findExpectedRank(const std::vector<int>& distr, const std::vector<int>& othersCumDistr)
{
    if (othersCumDistr.size() == 0) return 0;

    int othersMaxVal = othersCumDistr.size();

    int n = 0;
    int m = othersCumDistr.back();
    long long sum = 0;

    for (int i = 0; i < (int) distr.size(); ++i)
    {
        n += distr[i];
        sum += (long long) distr[i] * (i == 0 ? 0 : i > othersMaxVal ? m : othersCumDistr[i - 1]);
        sum -= (long long) distr[i] * (i >= othersMaxVal ? 0 : m - othersCumDistr[i]);
    }

    if (n == 0 || m == 0) return 0;

    long long combs = (long long) n * m;

    return (double) sum / combs;
}

void addToDistr(int val, std::vector<int>& distr)
{
    while (val >= (int) distr.size()) distr.push_back(0);
    ++distr[val];
}

bool isDone(const Occurs& occurs)
{
    return std::all_of(occurs.begin(), occurs.end(), [](int n) {return n == 0;});
}

int randNumRolls(Occurs& occurs, int extraDice, int goods)
{
    int rolls = 0;
    while (!isDone(occurs) && rolls < goods)
    {
        ++rolls;
        randRemOcc(extraDice, occurs);
    }
    return rolls;
}

int randNumNewHits(int left)
{
    int newHits = 0;
    for (int i = 0; i < left; ++i)
    {
        if (randDiceRoll() == 0) ++newHits;
    }
    return newHits;
}

// Optimized version of simColl for simulations
void miniSimColl(State& state, int num, int collected)
{
    bool cont = true;
    while (cont)
    {
        cont = false;
        int left = NUM_DICE - collected;
        Move mv = getCollMove(state.free, state.goods, num, left);

        switch (mv.id)
        {
        case M_STOP_COLL:
            break;
        
        case M_CONT_COLL:
            if (state.goods > 0 && left > 0)
            {
                cont = true;
                --state.goods;
                int newHits = randNumNewHits(left);
                collected += newHits;
            }
            else assert(false);
            break;

        default:
            assert(false);
        }
    }

    int reward = num * collected;
    state.score += reward;
}

// Optimized version of simMove for simulations
void miniSimMove(State& state, Occurs& occurs, int extraDice, int reward)
{
    int rolls = randNumRolls(occurs, extraDice, state.goods);
    bool won = isDone(occurs);

    state.goods -= rolls;
    if (won) state.score += reward;
}

// Optimized version of simTurn for simulations
void miniSimTurn(State& state)
{
    Occurs diceOccurs;

    state.goods += GOODS_PER_TURN;

    diceRoll:

    randOcc(diceOccurs);
    Move mv = getMove(state.free, state.goods, diceOccurs);

    int id = mv.id;
    switch (id)
    {
    case M_ERROR:
    case M_STOP_COLL:
    case M_CONT_COLL:
        assert(false);
        break;

    case M_REROLL:
        if (state.goods > 0 && !MISERE)
        {
            --state.goods;
            goto diceRoll;
        }
        else assert(false);
        break;
    
    default:
        if (id >= 0 && id < NUM_COMBOS && isFree(state.free, id))
        {
            state.free = setUsed(state.free, id);
            int num = combos[id]->getCollectNumber();
            if (num > 0)
            {
                int collected = diceOccurs[num - 1];
                miniSimColl(state, num, collected);
            }
            else
            {
                Target t = combos[id]->getTargetByArgs(mv.args);
                Occurs newOccurs = t.occurs;
                int extraDice = NUM_DICE - diceInOcc(newOccurs);
                remOcc(newOccurs, diceOccurs);
                miniSimMove(state, newOccurs, extraDice, t.points);
            }
        }
        else assert(false);
    }
}

int miniSimGame(State& state)
{
    while (state.free)
    {
        miniSimTurn(state);
    }
    return state.score + state.goods * POINTS_PER_GOOD;
}

int miniSimGameReroll(State& state)
{
    state.goods -= 1 + GOODS_PER_TURN;
    return miniSimGame(state);
}

int miniSimGameColl(State& state, int num, int collected)
{
    miniSimColl(state, num, collected);
    return miniSimGame(state);
}

int miniSimGameMove(State& state, const Occurs& occurs, int extraDice, int reward)
{
    Occurs newOccurs = occurs;
    miniSimMove(state, newOccurs, extraDice, reward);
    return miniSimGame(state);
}

int miniSimGameCollStop(State& state, int num, int collected)
{
    state.score += num * collected;
    return miniSimGame(state);
}

int miniSimGameCollCont(State& state, int num, int collected)
{
    --state.goods;
    int newHits = randNumNewHits(NUM_DICE - collected);
    collected += newHits;
    miniSimColl(state, num, collected);
    return miniSimGame(state);
}

void findOthersCumDistr(std::vector<State>& states, const Config& config)
{
    int numSims = config.competitive;
    std::vector<int>& othersCumDistr = states[config.player].othersCumDistr;
    othersCumDistr.clear();

    for (int other = 0; other < config.numPlayers; ++other)
    {
        if (other == config.player) continue;

        for (int i = 0; i < numSims; ++i)
        {
            State otherTempState = states[other];
            addToDistr(miniSimGame(otherTempState), othersCumDistr);
        }
    }

    for (int i = 1; i < (int) othersCumDistr.size(); ++i)
    {
        othersCumDistr[i] += othersCumDistr[i - 1];
    }
}

double getMean(const std::vector<int>& distr)
{
    long long cnt = 0;
    long long sum = 0;
    for (int i = 0; i < (int) distr.size(); ++i)
    {
        cnt += distr[i];
        sum += i * distr[i];
    }
    return (double) sum / cnt;
}

const double MAX_COMP_SLACK = 4;

Move getCompetitiveCollMove(const std::vector<State>& states, int num, int collected, const Config& config)
{
    int numSims = config.competitive;
    std::vector<int> distr;

    const State& state = states[config.player];
    int free = state.free;
    int goods = state.goods;

    Move best;
    double bestExpected = -1;
    for (Move& option : getCollMoveOptionsScored(free, goods, num, NUM_DICE - collected))
    {
        double expected = state.score + num * collected + option.score;

        if (bestExpected < 0) bestExpected = expected;
        else if (std::abs(bestExpected - expected) > MAX_COMP_SLACK) break;

        distr.clear();

        for (int i = 0; i < numSims; ++i)
        {
            int res = 0;
            State tempState = state;
            if (option.id == M_STOP_COLL)
            {
                res = miniSimGameCollStop(tempState, num, collected);
            }
            else if (option.id == M_CONT_COLL)
            {
                res = miniSimGameCollCont(tempState, num, collected);
            }
            addToDistr(res, distr);
        }

        option.score = findExpectedRank(distr, state.othersCumDistr);
        best = std::max(best, option);

        // std::cerr << " " << option.toString() << ": " << option.score << " / " << getMean(distr) << " / " << expected << std::endl;
    }

    return best;
}

Move getCompetitiveMove(const std::vector<State>& states, const Occurs& diceOccurs, const Config& config)
{
    int numSims = config.competitive;
    std::vector<int> distr;

    const State& state = states[config.player];
    int free = state.free;
    int goods = state.goods;

    Move best;
    double bestExpected = -1;
    for (Move& option : getMoveOptionsScored(free, goods, diceOccurs))
    {
        double expected = state.score + option.score;

        if (bestExpected < 0) bestExpected = expected;
        else if (std::abs(bestExpected - expected) > MAX_COMP_SLACK) break;

        int id = option.id;
        State newState = state;
        newState.free = id != M_REROLL ? setUsed(free, id) : free;
        int num = id != M_REROLL ? combos[id]->getCollectNumber() : 0;

        distr.clear();

        for (int i = 0; i < numSims; ++i)
        {
            int res = 0;
            State tempNewState = newState;
            if (id == M_REROLL)
            {
                res = miniSimGameReroll(tempNewState);
            }
            else if (num > 0)
            {
                int curr = diceOccurs[num - 1];
                res = miniSimGameColl(tempNewState, num, curr);
            }
            else
            {
                Target t = combos[id]->getTargetByArgs(option.args);
                Occurs newOccurs = t.occurs;
                int extraDice = NUM_DICE - diceInOcc(newOccurs);
                remOcc(newOccurs, diceOccurs);
                res = miniSimGameMove(tempNewState, newOccurs, extraDice, t.points);
            }
            addToDistr(res, distr);
        }

        option.score = findExpectedRank(distr, state.othersCumDistr);
        best = std::max(best, option);

        // std::cerr << " " << option.toString() << ": " << option.score << " / " << getMean(distr) << " / " << expected << std::endl;
    }

    return best;
}

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
    if (config.logMode != LOG_READ) std::cout << "Number of rolls to complete the combination (-1 for if you did not): ";
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

void simColl(std::vector<State>& states, int num, int collected, Config& config)
{
    State& state = states[config.player];

    bool first = true;
    bool cont = true;
    bool fail = false;
    while (cont)
    {
        cont = false;
        int left = NUM_DICE - collected;
        Move bestMv = getCollMove(state.free, state.goods, num, left);

        if (config.evalLM) state.updateExpScore(num * collected + bestMv.score, first ? R_NONE : R_LUCK);
        first = false;

        moveSelectColl:

        Move mv = config.manualMoves || config.logMode == LOG_READ ? chooseMove(config) :
                  (config.competitive ? getCompetitiveCollMove(states, num, collected, config) : bestMv);

        if (!config.manualMoves && config.verbose) std::cout << "Move: " << mv.toString() << std::endl;

        switch (mv.id)
        {
        case M_STOP_COLL:
            if (config.logMode == LOG_WRITE) config.logOut << mv.toString() << std::endl;

            if (config.evalLM) state.updateExpScore(num * collected + getScore(state.free, state.goods), R_MISTAKE);
            break;
        
        case M_CONT_COLL:
            if (state.goods > 0 && left > 0)
            {
                if (config.logMode == LOG_WRITE) config.logOut << mv.toString() << std::endl;

                if (config.evalLM) state.updateExpScore(num * collected + getCollContScore(state.free, state.goods, num, left), R_MISTAKE);
                cont = true;
                --state.goods;
                int newHits;
                if (config.manualRolls || config.logMode == LOG_READ) newHits = chooseNumNewHits(num, left, config);
                else newHits = randNumNewHits(left);

                if (config.logMode == LOG_WRITE) config.logOut << newHits << std::endl;

                if (!config.manualRolls && config.verbose)
                {
                    std::cout << "Number of " << numNames[num - 1] << " rolled: " << newHits << std::endl;
                }

                collected += newHits;
            }
            else fail = true;
            break;

        default:
            fail = true;
        }

        if (config.manualMoves && fail)
        {
            std::cout << "Invalid move." << std::endl;
            fail = false;
            goto moveSelectColl;
        }
    }

    assert(!fail);

    int reward = num * collected;
    state.score += reward;

    if (config.verbose)
    {
        std::cout << "Won " << reward << " points." << std::endl;
    }
}

void simMove(std::vector<State>& states, Occurs& occurs, int extraDice, int reward, Config& config)
{
    State& state = states[config.player];

    bool instaDone = isDone(occurs);
    int rolls;
    bool won;

    if (config.manualRolls || config.logMode == LOG_READ)
    {
        if (!instaDone && state.goods > 0) rolls = chooseNumRolls(config);
        else rolls = 0;

        won = instaDone || (rolls > 0 && rolls <= state.goods);
        if (!won) rolls = state.goods;
    }
    else
    {
        rolls = randNumRolls(occurs, extraDice, state.goods);
        won = isDone(occurs);
    }

    if (!instaDone && state.goods > 0 && config.logMode == LOG_WRITE) config.logOut << (won ? rolls : -1) << std::endl;

    if (!instaDone && !config.manualRolls && config.verbose)
    {
        if (won) std::cout << "Took " << rolls << " rolls to complete the combination. " << std::endl;
        else std::cout << "Did not complete the combination." << std::endl;
    }

    state.goods -= rolls;
    if (won) state.score += reward;

    if (config.evalLM) state.updateExpScore(getScore(state.free, state.goods), instaDone ? R_NONE : R_LUCK);

    if (won && config.verbose)
    {
        std::cout << "Won " << reward << " points." << std::endl;
    }
}

void simTurn(std::vector<State>& states, Config& config)
{
    State& state = states[config.player];
    Occurs diceOccurs;

    if (config.verbose)
    {
        std::cout << std::endl;
        if (config.numPlayers > 1) std::cout << "Player " << config.player + 1 << std::endl;
        std::cout << "Turn: " << state.getTurn() + 1 << "/" << NUM_COMBOS << std::endl;
        std::cout << "Current score: " << state.score << std::endl;
    }

    if (config.evalLM)
    {
        state.updateExpScore(getScore(state.free, state.goods), R_NONE);
        std::cout << "Expected final score: " << std::fixed << std::setprecision(IO_PRECISION) << state.expScore << std::endl;
    }

    if (!config.manualMoves && config.competitive) findOthersCumDistr(states, config);

    state.goods += GOODS_PER_TURN;

    diceRoll:

    if (config.manualRolls || config.logMode == LOG_READ) chooseOcc(diceOccurs, config);
    else randOcc(diceOccurs);

    if (config.logMode == LOG_WRITE) printOcc(diceOccurs, config, true);
    if (!config.manualRolls && config.verbose) printOcc(diceOccurs, config, false);

    Move bestMv = getMove(state.free, state.goods, diceOccurs);

    if (config.evalLM) state.updateExpScore(bestMv.score, R_LUCK);

    if (config.verbose)
    {
        std::cout << "Goods: " << state.goods << std::endl;
    }

    bool fail = false;

    moveSelect:

    Move mv = config.manualMoves || config.logMode == LOG_READ ? chooseMove(config) :
              (config.competitive ? getCompetitiveMove(states, diceOccurs, config) : bestMv);

    if (!config.manualMoves && config.verbose) std::cout << "Move: " << mv.toString() << std::endl;

    int id = mv.id;
    switch (id)
    {
    case M_ERROR:
    case M_STOP_COLL:
    case M_CONT_COLL:
        fail = true;
        break;

    case M_REROLL:
        if (state.goods > 0 && !MISERE)
        {
            if (config.logMode == LOG_WRITE) config.logOut << mv.toString() << std::endl;

            --state.goods;
            if (config.evalLM) state.updateExpScore(getContScore(state.free, state.goods), R_MISTAKE, bestMv.toString());
            goto diceRoll;
        }
        else fail = true;
        break;
    
    default:
        if (id >= 0 && id < NUM_COMBOS && isFree(state.free, id))
        {
            if (config.logMode == LOG_WRITE) config.logOut << mv.toString() << std::endl;

            state.free = setUsed(state.free, id);
            int num = combos[id]->getCollectNumber();
            if (num > 0)
            {
                int collected = diceOccurs[num - 1];
                if (config.evalLM) state.updateExpScore(num * collected + getCollScore(state.free, state.goods, num, NUM_DICE - collected), R_MISTAKE, bestMv.toString());
                simColl(states, num, collected, config);
            }
            else
            {
                Target t = combos[id]->getTargetByArgs(mv.args);
                Occurs newOccurs = t.occurs;
                int extraDice = NUM_DICE - diceInOcc(newOccurs);
                remOcc(newOccurs, diceOccurs);
                if (config.evalLM) state.updateExpScore(getMoveScore(state.free, state.goods, newOccurs, extraDice, t.points), R_MISTAKE, bestMv.toString());
                simMove(states, newOccurs, extraDice, t.points, config);
            }
        }
        else fail = true;
    }

    if (config.manualMoves && fail)
    {
        if (config.verbose) std::cout << "Invalid move" << std::endl;
        fail = false;
        goto moveSelect;
    }

    assert(!fail);
}

struct Result
{
    int score = 0;
    int rank = 0;
};

double displayRank(int rank, const Config& config)
{
    return (1 + config.numPlayers - rank) / 2.0;
}

std::vector<Result> simGame(Config& config)
{
    // For faster default simulations
    if (config.isDefault())
    {
        State state;
        int res = miniSimGame(state);
        return {{res, 0}};
    }

    if (config.verbose) std::cout << std::endl;

    std::vector<State> states(config.numPlayers);
    for (int turn = 0; turn < NUM_COMBOS; ++turn)
    {
        for (int player = 0; player < config.numPlayers; ++player)
        {
            config.setPlayer(player);
            simTurn(states, config);
        }
    }

    for (int player = 0; player < config.numPlayers; ++player)
    {
        states[player].score += states[player].goods * POINTS_PER_GOOD;
        states[player].goods = 0;

        if (config.evalLM) states[player].updateExpScore(0, R_NONE);
    }

    std::vector<Result> results(config.numPlayers);
    for (int player = 0; player < config.numPlayers; ++player)
    {
        results[player].score = states[player].score;
        for (int other = 0; other < config.numPlayers; ++other)
        {
            if (other == player) continue;
            if (states[other].score < states[player].score) ++results[player].rank;
            else if (states[other].score > states[player].score) --results[player].rank;
        }
    }

    for (int player = 0; player < config.numPlayers; ++player)
    {
        if (config.verbose) std::cout << std::endl;

        if (config.verbose)
        {
            if (config.numPlayers > 1) std::cout << "Player " << player + 1 << std::endl;
            std::cout << "Final score: " << results[player].score << std::endl;
            if (config.numPlayers > 1) std::cout << "Final rank: " << displayRank(results[player].rank, config) << std::endl;
        }

        if (config.evalLM)
        {
            states[player].printScoreByReason();
        }
    }

    return results;
}

struct Stats
{
    double mean = 0;
    double stdev = 0;
    double perc5 = 0;
    double perc25 = 0;
    double perc50 = 0;
    double perc75 = 0;
    double perc95 = 0;
    std::vector<double> modes;
    std::vector<int> distr;
};

const bool PRINT_RUN_MEAN = false;

Stats findStats(int n, Config& config, int povPlayer = 0, bool useRank = false)
{
    Stats stats;

    auto getVal = [&] (int r) {
        return !useRank ? r : displayRank(r - config.numPlayers, config);
    };

    double runSumScore = 0;
    double runSumRank = 0;

    for (int i = 0; i < n; ++i)
    {
        auto results = simGame(config);
        int r = !useRank ? results[povPlayer].score : results[povPlayer].rank + config.numPlayers;
        addToDistr(r, stats.distr);

        if (PRINT_RUN_MEAN)
        {
            runSumScore += results[povPlayer].score;
            runSumRank += displayRank(results[povPlayer].rank, config);
            std::cerr << i + 1 << ": " << runSumScore / (i + 1) << " " << runSumRank / (i + 1) << std::endl;
        }
    }

    double sum = 0;
    double sqSum = 0;
    int modeCnt = 0;
    int cumDistr = 0;

    for (int r = 0; r < (int) stats.distr.size(); ++r)
    {
        double val = getVal(r);

        sum += (long long) stats.distr[r] * val;
        sqSum += (long long) stats.distr[r] * val * val;

        if (stats.distr[r] > modeCnt)
        {
            stats.modes.clear();
            modeCnt = stats.distr[r];
        }

        if (stats.distr[r] == modeCnt) stats.modes.push_back(val);

        if (cumDistr < n / 20) stats.perc5 = val;
        if (cumDistr < n / 4) stats.perc25 = val;
        if (cumDistr < n / 2) stats.perc50 = val;
        if (cumDistr < 3 * n / 4) stats.perc75 = val;
        if (cumDistr < 19 * n / 20) stats.perc95 = val;

        cumDistr += stats.distr[r];
    }

    stats.mean = (double) sum / n;
    stats.stdev = sqrt((double) sqSum / n - stats.mean * stats.mean);

    return stats;
}

std::string modelName()
{
    return (!MISERE ? "normal_" : "misere_") + (NUM_TRIALS > 0 ? std::to_string(NUM_TRIALS) : "exact") + ".model";
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
            int numPlayers = 0;

            std::cout << "Number of players: ";
            while (numPlayers <= 0)
            {
                std::cin >> numPlayers;
            }

            std::vector<bool> manualRollsAll(numPlayers);
            std::vector<bool> manualMovesAll(numPlayers);
            std::vector<int> competitiveAll(numPlayers);
            bool evalLM;

            for (int player = 0; player < numPlayers; ++player)
            {
                int temp;

                if (numPlayers > 1) std::cout << "Player " << player + 1 << ": ";
                std::cout << "Manual rolls (0 / 1): ";
                std::cin >> temp;
                manualRollsAll[player] = temp;

                if (numPlayers > 1) std::cout << "Player " << player + 1 << ": ";
                std::cout << "Manual moves (0 / 1): ";
                std::cin >> temp;
                manualMovesAll[player] = temp;

                if (numPlayers > 1 && !manualMovesAll[player])
                {
                    if (numPlayers > 1) std::cout << "Player " << player + 1 << ": ";
                    std::cout << "Competitive (0 / 1+): ";
                    std::cin >> competitiveAll[player];
                }
                else competitiveAll[player] = 0;
            }

            std::cout << "Evaluate luck and mistakes (0 / 1): ";
            std::cin >> evalLM;

            lastLogName = std::to_string(timestamp()) + ".log";

            Config config(numPlayers, true, evalLM, LOG_WRITE, lastLogName, manualRollsAll, manualMovesAll, competitiveAll);
            simGame(config);
        }
        else if (cmd == "replay")
        {
            bool evalLM;
            std::string logName;

            std::cout << (lastLogName != "" ? "Log name or 'last' for last one: " : "Log name: ");
            std::cin >> logName;

            if (lastLogName != "" && logName == "last") logName = lastLogName;

            std::cout << "Evaluate luck and mistakes (0 / 1): ";
            std::cin >> evalLM;

            Config config(0, true, evalLM, LOG_READ, logName);
            simGame(config);
        }
        else if (cmd == "test")
        {
            int numPlayers = 0;
            int povPlayer = -1;
            bool useRank;

            std::cout << "Number of players: ";
            while (numPlayers <= 0)
            {
                std::cin >> numPlayers;
            }

            std::vector<int> competitiveAll;

            if (numPlayers == 1)
            {
                povPlayer = 0;
                useRank = false;
            }
            else
            {
                std::cout << "PoV player: ";
                while (povPlayer < 0 || povPlayer >= numPlayers)
                {
                    std::cin >> povPlayer;
                    --povPlayer;
                }

                std::cout << "Score or rank (0 / 1): ";
                std::cin >> useRank;

                competitiveAll.resize(numPlayers);

                for (int player = 0; player < numPlayers; ++player)
                {
                    if (numPlayers > 1) std::cout << "Player " << player + 1 << ": ";
                    std::cout << "Competitive (0 / 1+): ";
                    std::cin >> competitiveAll[player];
                }
            }

            Config config(numPlayers, false, false, LOG_NONE, "", {}, {}, competitiveAll);

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

            Stats stats = findStats(numTests, config, povPlayer, useRank);

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
            std::cout << std::endl;
            std::cout << "Expected score: " << getScore(NUM_MASKS - 1, 0) << std::endl;
        }
        if (cmd == "credits")
        {
            std::cout << std::endl;
            std::cout << credits << std::endl;
        }
        if (cmd == "help")
        {
            std::cout << std::endl;
            std::cout << help << std::endl;
        }
        if (cmd == "exit") break;
    }
}

int main()
{
    generator.seed(time(nullptr));

    std::cout << "Normal play or misere play (0 / 1): ";
    std::cin >> MISERE;

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
