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

const double EPS = 1e-9;
const double INF = 1e9;

// Game configuration
const int NUM_DICE = 6;
const int NUM_SIDES = 6;
const int NUM_COMBOS = 14;
const int GOODS_PER_TURN = 5;
const int POINTS_PER_GOOD = 1;

// Depends on game configuration
const int ZERO_CODE = 0;
const int ZERO_ORD_CODE = 0;
const int VALID_CODES = 30; // Sum_{0 <= k <= NUM_DICE} #ways to partition k into at most NUM_SIDES partitions
const int STARTS_ORD_CODES[NUM_DICE + 2] = {0, 1, 7, 28, 84, 210, 462, 924}; // 0 <= t <= NUM_DICE + 1: Sum_{0 <= k <= t} #ways to partition k into NUM_SIDES ordered partitions
const int VALID_ORD_CODES = STARTS_ORD_CODES[NUM_DICE + 1];
const int START_FULL_ORD_CODES = STARTS_ORD_CODES[NUM_DICE];
const int NUM_MASKS = 1 << NUM_COMBOS;
const int INITIAL_MASK = NUM_MASKS - 1;
const int MAX_GOODS = NUM_COMBOS * GOODS_PER_TURN;

int NUM_TRIALS; // Trials per state
bool MISERE; // Normal play or misere play

double worstScore()
{
    return !MISERE ? -INF : INF;
}

template <class T>
bool smartLt(T a, T b)
{
    return !MISERE ? a < b : b < a;
}

using Occurs = std::array<int, NUM_SIDES>;

std::unordered_map<int, int> rawCodeToCode;
std::unordered_map<int, int> ordRawCodeToOrdCode;

int occursToCode(Occurs occurs, bool ordered)
{
    std::unordered_map<int, int>& currRawCodeToCode = !ordered ? rawCodeToCode : ordRawCodeToOrdCode;
    if (!ordered) std::sort(occurs.begin(), occurs.end(), std::greater<int>());

    int rawCode = 0;
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        rawCode = rawCode * (NUM_DICE + 1) + occurs[i];
    }
    return currRawCodeToCode.insert({rawCode, currRawCodeToCode.size()}).first->second;
}

int fact(int n)
{
    return n ? n * fact(n - 1) : 1;
}

int diceInOccurs(const Occurs& occurs)
{
    return std::accumulate(occurs.begin(), occurs.end(), 0);
}

double occursProb(const Occurs& occurs)
{
    int numDice = diceInOccurs(occurs);
    double prob = fact(numDice);
    for (int cnt : occurs)
    {
        prob /= fact(cnt);
    }
    return prob / pow(NUM_SIDES, numDice);
}

int diceInCode[VALID_CODES];
int diceInOrdCode[VALID_ORD_CODES];
Occurs codeOccurs[VALID_CODES];
Occurs ordCodeOccurs[VALID_ORD_CODES];
double ordCodeProb[VALID_ORD_CODES];

void genCodesRec(Occurs& occurs, int pos, int dice)
{
    if (pos == NUM_SIDES && dice == 0)
    {
        int code = occursToCode(occurs, false);
        int ordCode = occursToCode(occurs, true);

        diceInCode[code] = diceInOccurs(occurs);
        diceInOrdCode[ordCode] = diceInCode[code];
        codeOccurs[code] = occurs;
        ordCodeOccurs[ordCode] = occurs;
        ordCodeProb[ordCode] = occursProb(occurs);
    }

    if (pos == NUM_SIDES) return;

    for (int i = 0; i <= dice; ++i)
    {
        occurs[pos] = i;
        genCodesRec(occurs, pos + 1, dice - i);
    }
}

Occurs remOccurs(const Occurs& occurs, const Occurs& diceOccurs)
{
    Occurs res;
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        res[i] = std::max(occurs[i] - diceOccurs[i], 0);
    }
    return res;
}

// Both result in code
int codeRemOrdCode[VALID_CODES][VALID_ORD_CODES];
int ordCodeRemOrdCode[VALID_ORD_CODES][VALID_ORD_CODES];

void genCodes()
{
    Occurs occurs;
    for (int i = 0; i <= NUM_DICE; ++i)
    {
        genCodesRec(occurs, 0, i);
    }

    for (int code = 0; code < VALID_CODES; ++code)
    {
        for (int diceOrdCode = 0; diceOrdCode < VALID_ORD_CODES; ++diceOrdCode)
        {
            codeRemOrdCode[code][diceOrdCode] = occursToCode(remOccurs(codeOccurs[code], ordCodeOccurs[diceOrdCode]), false);
        }
    }

    for (int ordCode = 0; ordCode < VALID_ORD_CODES; ++ordCode)
    {
        for (int diceOrdCode = 0; diceOrdCode < VALID_ORD_CODES; ++diceOrdCode)
        {
            ordCodeRemOrdCode[ordCode][diceOrdCode] = occursToCode(remOccurs(ordCodeOccurs[ordCode], ordCodeOccurs[diceOrdCode]), false);
        }
    }
}

#define NONUM -1

struct Target
{
    int ordCode;
    int points;
    std::vector<int> args;

    Target(int ordCode, int points, const std::vector<int>& args):
        ordCode(ordCode),
        points(points),
        args(args) {}

    Target(int ordCode, int points):
        Target(ordCode, points, {}) {}
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
        assert(false);
        return 0;
    }

    virtual int getCollectNumber() const
    {
        return NONUM;
    }

    virtual const std::vector<Target>& getTargets(int diceOrdCode) const
    {
        assert(false);
        static std::vector<Target> targets;
        return targets;
    }

    virtual Target getTargetByArgs(const std::vector<int>& args) const
    {
        assert(false);
        return Target(ZERO_ORD_CODE, 0);
    }

    virtual ~Combo() {}

protected:

    const std::string name;
    const int points;
};

struct CollectCombo : Combo
{
    CollectCombo(const std::string& name, int number):
        Combo(name, 0),
        number(number) {}

    int getCollectNumber() const override
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
        targets({Target(occursToCode(occurs, true), points)}) {}

    int getNumArgs() const override
    {
        return 0;
    }

    const std::vector<Target>& getTargets(int diceOrdCode) const override
    {
        return targets;
    }

    Target getTargetByArgs(const std::vector<int>& args) const override
    {
        return targets.front();
    }

protected:

    std::vector<Target> targets;
};

bool checkArgs(const std::vector<int>& args)
{
    std::unordered_set<int> takenArgs;
    for (int arg : args)
    {
        if (arg < 0 || arg >= NUM_SIDES) return false;
        auto res = takenArgs.insert(arg);
        if (!res.second) return false;
    }
    return true;
}

struct PermCombo : Combo
{
    PermCombo(const std::string& name, int points, const Occurs& occurs):
        PermCombo(name, points, occurs, false) {}

    int getNumArgs() const override
    {
        return numArgs;
    }

    const std::vector<Target>& getTargets(int diceOrdCode) const override
    {
        return cache[diceOrdCode];
    }

    Target getTargetByArgs(const std::vector<int>& args) const override;

protected:

    int templateCode;

    int numArgs;
    std::vector<Target> cache[VALID_ORD_CODES];

    PermCombo(const std::string& name, int points, const Occurs& occurs, bool isChild):
        Combo(name, points),
        templateCode(occursToCode(occurs, false))
    {
        for (numArgs = 0; numArgs < NUM_SIDES && codeOccurs[templateCode][numArgs] > 0; ++numArgs);

        if (!isChild) buildCache();
    }

    void buildCache()
    {
        for (int diceOrdCode = START_FULL_ORD_CODES; diceOrdCode < VALID_ORD_CODES; ++diceOrdCode)
        {
            cache[diceOrdCode] = rawGetTargets(diceOrdCode);
        }
    }

    virtual int getPoints(int ordCode) const
    {
        return points;
    }

    int makeOrdCodeByArgs(const std::vector<int>& args) const
    {
        Occurs occurs;
        occurs.fill(0);
        for (int i = 0; i < numArgs; ++i)
        {
            occurs[args[i]] = codeOccurs[templateCode][i];
        }
        return occursToCode(occurs, true);
    }

    std::vector<Target> rawGetTargets(int diceOrdCode) const
    {
        std::vector<std::pair<Target, int>> possTargetCodes; 

        std::vector<int> args(numArgs, 0);
        while (true)
        {
            if (checkArgs(args))
            {
                Target target = getTargetByArgs(args);
                int code = ordCodeRemOrdCode[target.ordCode][diceOrdCode];
                possTargetCodes.push_back({target, code});
            }

            bool inc = true;
            int idx = numArgs - 1;
            while (inc && idx >= 0)
            {
                ++args[idx];
                if (args[idx] == NUM_SIDES)
                {
                    args[idx] = 0;
                    inc = true;
                    --idx;
                }
                else inc = false;
            }
            if (inc) break;
        }

        std::vector<Target> targets;

        for (int i = 0; i < (int) possTargetCodes.size(); ++i)
        {
            auto& [target, code] = possTargetCodes[i];

            bool keep = true;
            for (int j = 0; j < (int) possTargetCodes.size() && keep; ++j)
            {
                if (j == i) continue;

                auto& [otherTarget, otherCode] = possTargetCodes[j];

                int cmp = 0;

                if (smartLt(otherTarget.points, target.points)) continue;

                for (int k = 0; k < NUM_SIDES && cmp != 1; ++k)
                {
                    if (codeOccurs[code][k] == codeOccurs[otherCode][k]) continue;
                    cmp = smartLt(codeOccurs[code][k], codeOccurs[otherCode][k]) ? 1 : -1;
                }

                if (cmp == 1) continue;
                else if (cmp == -1) keep = false;
                else if (cmp == 0)
                {
                    if (smartLt(target.points, otherTarget.points)) keep = false;
                    else keep = i < j;
                }
            }

            if (keep) targets.push_back(target);
        }

        return targets;
    }
};

struct SPermCombo : public PermCombo
{
    SPermCombo(const std::string& name, const Occurs& occurs):
        PermCombo(name, 0, occurs, true)
    {
        buildCache();
    }

protected:

    int getPoints(int ordCode) const override
    {
        int points = 0;
        for (int i = 0; i < NUM_SIDES; ++i)
        {
            points += (i + 1) * ordCodeOccurs[ordCode][i];
        }
        return points;
    }
};

Target PermCombo::getTargetByArgs(const std::vector<int>& args) const
{
    int ordCode = makeOrdCodeByArgs(args);
    int currPoints = getPoints(ordCode);
    return Target(ordCode, currPoints, args);
}

// Game configuration
std::array<Combo*, NUM_COMBOS> combos;

void setCombos()
{
    combos = {
        new CollectCombo("Collect", 0),
        new CollectCombo("Collect", 1),
        new CollectCombo("Collect", 2),
        new CollectCombo("Collect", 3),
        new CollectCombo("Collect", 4),
        new CollectCombo("Collect", 5),
        new   SPermCombo("Three pairs",        {2, 2, 2, 0, 0, 0}),
        new   SPermCombo("Two triples",        {3, 3, 0, 0, 0, 0}),
        new    PermCombo("Four of a kind", 40, {4, 0, 0, 0, 0, 0}),
        new   FixedCombo("Kamerun",        45, {0, 0, 0, 1, 2, 3}),
        new   FixedCombo("Straight",       50, {1, 1, 1, 1, 1, 1}),
        new    PermCombo("Six of a kind",  60, {6, 0, 0, 0, 0, 0}),
        new   FixedCombo("General",        70, {0, 0, 0, 0, 0, 6}),
        new   FixedCombo("Despirala",      80, {5, 0, 0, 0, 0, 1})
    };
}

const int MAX_EXTRA_DICE = 2; // NUM_DICE - min dice in combo

std::mt19937 generator;
std::uniform_int_distribution<int> diceDistribution(0, NUM_SIDES - 1);

int randDiceRoll()
{
    return diceDistribution(generator);
}

Occurs randOccurs(int numDice)
{
    Occurs occurs;
    occurs.fill(0);
    for (int i = 0; i < numDice; ++i)
    {
        ++occurs[randDiceRoll()];
    }
    return occurs;
}

int randOrdCode(int numDice = NUM_DICE)
{
    return occursToCode(randOccurs(numDice), true);
}

double rollsDistr[MAX_EXTRA_DICE + 1][VALID_CODES][MAX_GOODS + 1];

void findRollsDistrSingle(int extraDice, int code)
{
    double changeProb = 0;
    int numDice = extraDice + diceInCode[code];
    for (int diceOrdCode = STARTS_ORD_CODES[numDice]; diceOrdCode < STARTS_ORD_CODES[numDice + 1]; ++diceOrdCode)
    {
        int newCode = codeRemOrdCode[code][diceOrdCode];
        if (newCode == code) continue;
        changeProb += ordCodeProb[diceOrdCode];
        for (int j = 0; j < MAX_GOODS; ++j)
        {
            rollsDistr[extraDice][code][j] += rollsDistr[extraDice][newCode][j] * ordCodeProb[diceOrdCode];
        }
    }
    for (int i = MAX_GOODS; i >= 0; --i)
    {
        rollsDistr[extraDice][code][i] = 0;
        for (int j = i - 1; j >= 0 ; --j)
        {
            rollsDistr[extraDice][code][i] += rollsDistr[extraDice][code][j] * pow(1 - changeProb, i - j - 1);
        }
    }
}

void findRollsDistr()
{
    for (int extraDice = 0; extraDice <= MAX_EXTRA_DICE; ++extraDice)
    {
        rollsDistr[extraDice][0][0] = 1;
        for (int code = 1; code < VALID_CODES; ++code)
        {
            if (extraDice + diceInCode[code] > NUM_DICE) continue; 
            findRollsDistrSingle(extraDice, code);
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
#define M_LIST_OPTIONS -99
#define M_ERROR -100

struct Move
{
    int id;
    Target target;
    double score;

    Move(int id, Target target, double score):
        id(id),
        target(target),
        score(score) {}

    Move(int id, double score):
        Move(id, Target(ZERO_ORD_CODE, 0), score) {}

    Move(int id, Target target):
        Move(id, target, worstScore()) {}

    Move(int id):
        Move(id, worstScore()) {}

    Move():
        Move(M_ERROR) {}

    Move(const std::string& name, const std::vector<int>& args);
    
    std::string toString(bool printTemplate = false) const;
};

bool operator<(const Move& a, const Move& b)
{
    return smartLt(a.score, b.score);
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

std::string Move::toString(bool printTemplate) const
{
    std::string s;
    switch (id)
    {
    case M_ERROR:
        s = "Error";
        break;

    case M_LIST_OPTIONS:
        s = "Options";
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
        int num = combos[id]->getCollectNumber();
        if (num != NONUM) s += " " + std::to_string(num + 1);
        else
        {
            for (int i = 0; i < (int) target.args.size(); ++i)
            {
                s += " ";
                if (!printTemplate) s += std::to_string(target.args[i] + 1);
                else s += 'A' + i;
            }
        }
    }
    return s;
}

Move::Move(const std::string& name, const std::vector<int>& args):
    Move()
{
    if (!checkArgs(args)) return;

    int nArgs = args.size();

    if ((name == "list options" || name == "options" || name == "list") && nArgs == 0) id = M_LIST_OPTIONS;
    else if ((name == "stop collecting" || name == "stop") && nArgs == 0) id = M_STOP_COLL;
    else if ((name == "continue collecting" || name == "continue") && nArgs == 0) id = M_CONT_COLL;
    else if (name == "reroll" && nArgs == 0) id = M_REROLL;
    else
    {
        for (int i = 0; i < NUM_COMBOS; ++i)
        {
            int num = combos[i]->getCollectNumber();
            if (name == stringToLower(combos[i]->getName()) && ((num != NONUM && nArgs == 1) || (num == NONUM && nArgs == combos[i]->getNumArgs())))
            {
                if (num == NONUM || args.front() == num)
                {
                    id = i;
                    if (num == NONUM) target = combos[i]->getTargetByArgs(args);
                    break;
                }
            }
        }
    }
}

bool isFree(int free, int idx)
{
    return (free >> idx) & 1;
}

int setUsed(int free, int idx)
{
    return free - (1 << idx);
}

std::vector<Move> getMoveOptions(int free, int goods, int diceOrdCode)
{
    std::vector<Move> options;
    for (int i = 0; i < NUM_COMBOS; ++i)
    {
        if (!isFree(free, i)) continue;
        int num = combos[i]->getCollectNumber();
        if (num != NONUM)
        {
            options.emplace_back(i);
        }
        else
        {
            for (const Target& target : combos[i]->getTargets(diceOrdCode))
            {
                options.emplace_back(i, target);
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
        sc += leftDistr[left][i] * ((left - i) * (num + 1) + getCollScore(free, goods - 1, num, i));
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

double getMoveScore(int free, int goods, int code, int extraDice, int points)
{
    double score = 0;
    double succP = 0;
    for (int i = 0; i <= goods; ++i)
    {
        double currP = rollsDistr[extraDice][code][i];
        if (currP > 0)
        {
            succP += currP;
            score += currP * getScore(free, goods - i);
        }
    }
    score += succP * points;
    score += (1 - succP) * getScore(free, 0);
    return score;
}

std::vector<Move> getMoveOptionsScored(int free, int goods, int diceOrdCode)
{
    std::vector<Move> options = getMoveOptions(free, goods, diceOrdCode);
    for (Move& option : options)
    {
        int id = option.id;
        int newFree = id != M_REROLL ? setUsed(free, id) : free;
        int num = id != M_REROLL ? combos[id]->getCollectNumber() : 0;

        if (id == M_REROLL)
        {
            option.score = getContScore(newFree, goods - 1);
        }
        else if (num != NONUM)
        {
            int collected = ordCodeOccurs[diceOrdCode][num];
            option.score = collected * (num + 1) + getCollScore(newFree, goods, num, NUM_DICE - collected);
        }
        else
        {
            int newCode = ordCodeRemOrdCode[option.target.ordCode][diceOrdCode];
            int extraDice = NUM_DICE - diceInOrdCode[option.target.ordCode];
            option.score = getMoveScore(newFree, goods, newCode, extraDice, option.target.points);
        }
    }
    std::sort(options.rbegin(), options.rend());
    return options;
}

// Optimized version of getMoveOptionsScored(...)[0]
Move getMove(int free, int goods, int diceOrdCode)
{
    return getMoveOptionsScored(free, goods, diceOrdCode)[0];

    Move best;
    for (int i = 0; i < NUM_COMBOS; ++i)
    {
        if (!isFree(free, i)) continue;
        int newFree = setUsed(free, i);
        int num = combos[i]->getCollectNumber();
        if (num != NONUM)
        {
            int collected = ordCodeOccurs[diceOrdCode][num];
            Move option(i, collected * (num + 1) + getCollScore(newFree, goods, num, NUM_DICE - collected));
            best = std::max(best, option);
        }
        else
        {
            for (const Target& target : combos[i]->getTargets(diceOrdCode))
            {
                int newCode = ordCodeRemOrdCode[target.ordCode][diceOrdCode];
                int extraDice = NUM_DICE - diceInOrdCode[target.ordCode];
                Move option(i, target, getMoveScore(newFree, goods, newCode, extraDice, target.points));
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
            Move move = getMove(free, goods, randOrdCode());
            score[free][goods] += move.score / NUM_TRIALS;
        }
    }
    else
    {
        for (int diceOrdCode = START_FULL_ORD_CODES; diceOrdCode < VALID_ORD_CODES; ++diceOrdCode)
        {
            Move move = getMove(free, goods, diceOrdCode);
            score[free][goods] += move.score * ordCodeProb[diceOrdCode];
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
    return getScore(INITIAL_MASK, 0);
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

    std::vector<double> othersCumDistr;

    State():
        free(INITIAL_MASK),
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
            std::cout.setf(std::ios::showpos);
            if (reason == R_LUCK && mult * delta > 0) std::cout << "Good luck: ";
            else if (reason == R_LUCK && mult * delta < 0) std::cout << "Bad luck: ";
            else if (reason == R_MISTAKE && -mult * delta < 1) std::cout << "Inaccuracy: ";
            else if (reason == R_MISTAKE && -mult * delta < 4) std::cout << "Mistake: ";
            else if (reason == R_MISTAKE && -mult * delta < 10) std::cout << "Blunder: ";
            else if (reason == R_MISTAKE && -mult * delta >= 10) std::cout << "Massive blunder: ";
            else if (reason == R_NONE) std::cout << "(Warning) No reason: ";
            else std::cout << "(Error) Invalid reason " << reason << ": ";
            std::cout << delta << std::endl;
            if (reason == R_MISTAKE && bestMvName != "") std::cout << "Best move was: " << bestMvName << std::endl;
            std::cout.unsetf(std::ios::showpos);
        }
        scoreByReason[reason] += delta;
        expScore = newExpScore;
    }

    void printScoreByReason()
    {
        std::cout << "Baseline score: " << getInitialScore() << std::endl;

        std::cout.setf(std::ios::showpos);
        std::cout << "Score due to luck: " << scoreByReason[R_LUCK] << std::endl;
        std::cout.unsetf(std::ios::showpos);

        std::cout << "Score due to mistakes: " << scoreByReason[R_MISTAKE] << std::endl;

        if (fabs(scoreByReason[R_NONE]) > IO_EPS)
        {
            std::cout.setf(std::ios::showpos);
            std::cout << "(Warning) Score for no reason: " << std::showpos << scoreByReason[R_NONE] << std::endl;
            std::cout.unsetf(std::ios::showpos);
        }
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

    template <class T>
    std::vector<T> vectorOrDefault(std::vector<T> vec)
    {
        assert(vec.empty() || (int) vec.size() == numPlayers);

        if (vec.empty()) return std::vector<T>(numPlayers, T());
        else return vec;
    }
};

double findExpectedRank(const std::vector<double>& distr, const std::vector<double>& othersCumDistr)
{
    double expRank = 0;
    double m = othersCumDistr.back();
    int othersMaxVal = othersCumDistr.size();
    for (int i = 0; i < (int) distr.size(); ++i)
    {
        expRank += distr[i] * (i == 0 ? 0 : i > othersMaxVal ? m : othersCumDistr[i - 1]);
        expRank -= distr[i] * (i >= othersMaxVal ? 0 : m - othersCumDistr[i]);
    }
    return expRank;
}

void addToDistr(int val, std::vector<int>& distr)
{
    while (val >= (int) distr.size()) distr.push_back(0);
    ++distr[val];
}

bool isDone(int code)
{
    return code == ZERO_CODE;
}

int randNumRolls(int& code, int extraDice, int goods)
{
    int rolls = 0;
    while (!isDone(code) && rolls < goods)
    {
        ++rolls;
        code = codeRemOrdCode[code][randOrdCode(diceInCode[code] + extraDice)];
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
        Move move = getCollMove(state.free, state.goods, num, left);

        switch (move.id)
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

    int points = collected * (num + 1);
    state.score += points;
}

// Optimized version of simMove for simulations
void miniSimMove(State& state, int code, int extraDice, int points)
{
    int rolls = randNumRolls(code, extraDice, state.goods);
    bool won = isDone(code);

    state.goods -= rolls;
    if (won) state.score += points;
}

// Optimized version of simTurn for simulations
void miniSimTurn(State& state)
{
    int diceOrdCode;

    state.goods += GOODS_PER_TURN;

    diceRoll:

    diceOrdCode = randOrdCode();
    Move move = getMove(state.free, state.goods, diceOrdCode);

    int id = move.id;
    switch (id)
    {
    case M_ERROR:
    case M_LIST_OPTIONS:
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
            if (num != NONUM)
            {
                int collected = ordCodeOccurs[diceOrdCode][num];
                miniSimColl(state, num, collected);
            }
            else
            {
                int newCode = ordCodeRemOrdCode[move.target.ordCode][diceOrdCode];
                int extraDice = NUM_DICE - diceInOrdCode[move.target.ordCode];
                miniSimMove(state, newCode, extraDice, move.target.points);
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

int miniSimGameMove(State& state, int code, int extraDice, int points)
{
    miniSimMove(state, code, extraDice, points);
    return miniSimGame(state);
}

int miniSimGameCollStop(State& state, int num, int collected)
{
    state.score += collected * (num + 1);
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

std::vector<double> fixMean(const std::vector<int>& q, double mu)
{
    int minVal = -1;
    int maxVal = q.size();
    int n = 0;
    for (int k = 0; k < maxVal; ++k)
    {
        if (minVal == -1 && q[k] > 0) minVal = k;
        n += q[k];
    }

    int stepsLeft = 20;

    double lamb1 = n / 2;
    double lamb2 = 0;
    double f;
    double g;

    if (mu - EPS <= minVal || mu + EPS >= maxVal - 1)
    {
        lamb1 = n;
        lamb2 = 0;
        goto exit;
    }

    do
    {
        f = 1;
        g = mu + 1;
        double jf1 = 0;
        double jf2 = 0;
        double jg1 = 0;
        double jg2 = 0;
        for (int k = minVal; k < maxVal; ++k)
        {
            int t = k + 1;
            double rec = 1 / (lamb1 + t * lamb2);
            f -= q[k] * rec;
            g -= t * q[k] * rec;
            jf1 += q[k] * rec * rec;
            jf2 += t * q[k] * rec * rec;
            jg1 += t * q[k] * rec * rec;
            jg2 += t * t * q[k] * rec * rec;
        }
        double jdet = jf1 * jg2 - jf2 * jg1;
        double ji1f = jg2 / jdet;
        double ji1g = -jf2 / jdet;
        double ji2f = -jg1 / jdet;
        double ji2g = jf1 / jdet;
        double d1 = f * ji1f + g * ji1g;
        double d2 = f * ji2f + g * ji2g;
        lamb1 -= d1;
        lamb2 -= d2;
    }
    while ((std::abs(f) > EPS || std::abs(g) > EPS) && stepsLeft-- > 0);

    if (!(std::abs(f) <= EPS && std::abs(g) <= EPS))
    {
        lamb1 = n;
        lamb2 = 0;
    }

    exit:

    std::vector<double> p(maxVal, 0);
    for (int k = minVal; k < maxVal; ++k)
    {
        int t = k + 1;
        p[k] = q[k] / (lamb1 + t * lamb2);
    }
    return p;
}

template <class T>
double getMean(const std::vector<T>& distr)
{
    double cnt = 0;
    double sum = 0;
    for (int i = 0; i < (int) distr.size(); ++i)
    {
        cnt += distr[i];
        sum += i * distr[i];
    }
    return sum / cnt;
}

void findOthersCumDistr(std::vector<State>& states, const Config& config)
{
    int numSims = config.competitive;
    std::vector<double>& othersCumDistr = states[config.player].othersCumDistr;
    othersCumDistr.clear();

    for (int other = 0; other < config.numPlayers; ++other)
    {
        if (other == config.player) continue;

        State& otherState = states[other];

        std::vector<int> empDistr;
        double expected = otherState.score + getScore(otherState.free, otherState.goods);
        for (int i = 0; i < numSims; ++i)
        {
            State otherTempState = otherState;
            addToDistr(miniSimGame(otherTempState), empDistr);
        }
        std::vector<double> distr = fixMean(empDistr, expected);

        // std::cerr << " Other: " << getMean(empDistr) << " / " << getMean(distr) << " / " << expected << std::endl;

        othersCumDistr.resize(std::max(othersCumDistr.size(), distr.size()), 0);
        double runSum = 0;
        for (int i = 0; i < (int) distr.size(); ++i)
        {
            runSum += distr[i];
            othersCumDistr[i] += runSum;
        }
    }
}

const double MAX_COMP_SLACK = 4;

Move getCompetitiveCollMove(const std::vector<State>& states, int num, int collected, const Config& config)
{
    int numSims = config.competitive;

    const State& state = states[config.player];
    int free = state.free;
    int goods = state.goods;

    Move best;
    double bestExpected = -1;
    for (Move& option : getCollMoveOptionsScored(free, goods, num, NUM_DICE - collected))
    {
        double expected = state.score + collected * (num + 1) + option.score;

        if (bestExpected < 0) bestExpected = expected;
        else if (std::abs(bestExpected - expected) > MAX_COMP_SLACK) break;

        std::vector<int> empDistr;
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
            addToDistr(res, empDistr);
        }

        std::vector<double> distr = fixMean(empDistr, expected);

        option.score = findExpectedRank(distr, state.othersCumDistr);
        best = std::max(best, option);

        // std::cerr << " " << option.toString() << ": " << option.score << " / " << getMean(empDistr) << " / " << getMean(distr) << " / " << expected << std::endl;
    }

    return best;
}

Move getCompetitiveMove(const std::vector<State>& states, int diceOrdCode, const Config& config)
{
    int numSims = config.competitive;

    const State& state = states[config.player];
    int free = state.free;
    int goods = state.goods;

    Move best;
    double bestExpected = -1;
    for (Move& option : getMoveOptionsScored(free, goods, diceOrdCode))
    {
        double expected = state.score + option.score;

        if (bestExpected < 0) bestExpected = expected;
        else if (std::abs(bestExpected - expected) > MAX_COMP_SLACK) break;

        int id = option.id;
        State newState = state;
        newState.free = id != M_REROLL ? setUsed(free, id) : free;
        int num = id != M_REROLL ? combos[id]->getCollectNumber() : 0;

        std::vector<int> empDistr;
        for (int i = 0; i < numSims; ++i)
        {
            int res = 0;
            State tempNewState = newState;
            if (id == M_REROLL)
            {
                res = miniSimGameReroll(tempNewState);
            }
            else if (num != NONUM)
            {
                int collected = ordCodeOccurs[diceOrdCode][num];
                res = miniSimGameColl(tempNewState, num, collected);
            }
            else
            {
                int newCode = ordCodeRemOrdCode[option.target.ordCode][diceOrdCode];
                int extraDice = NUM_DICE - diceInOrdCode[option.target.ordCode];
                res = miniSimGameMove(tempNewState, newCode, extraDice, option.target.points);
            }
            addToDistr(res, empDistr);
        }

        std::vector<double> distr = fixMean(empDistr, expected);

        option.score = findExpectedRank(distr, state.othersCumDistr);
        best = std::max(best, option);

        // std::cerr << " " << option.toString() << ": " << option.score << " / " << getMean(empDistr) << " / " << getMean(distr) << " / " << expected << std::endl;
    }

    return best;
}

void printOrdCode(int ordCode, Config& config, bool useConfig)
{
    if (!useConfig || config.logMode != LOG_WRITE) std::cout << "Rolled: ";
    std::ostream& out = useConfig && config.logMode == LOG_WRITE ? config.logOut : std::cout;

    bool first = true;
    for (int i = 0; i < NUM_SIDES; ++i)
    {
        for (int j = 0; j < ordCodeOccurs[ordCode][i]; ++j)
        {
            if (!first) out << " ";
            out<< i + 1;
            first = false;
        }
    }
    out << std::endl;
}

int chooseOrdCode(Config& config)
{
    if (config.logMode != LOG_READ) std::cout << "Rolled: ";
    std::istream& in = config.logMode == LOG_READ ? config.logIn : std::cin;

    Occurs occurs;
    occurs.fill(0);
    int side;
    for (int i = 0; i < NUM_DICE; ++i)
    {
        in >> side;
        if (side >= 1 && side <= NUM_SIDES) ++occurs[side - 1];
        else --i;
    }

    return occursToCode(occurs, true);
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
        if (isNumber(word)) args.push_back(stoi(word) - 1);
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
    if (config.logMode != LOG_READ) std::cout << "Number of " << numNames[num] << " rolled: ";
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
        Move bestMove = getCollMove(state.free, state.goods, num, left);

        if (config.evalLM) state.updateExpScore(collected * (num + 1) + bestMove.score, first ? R_NONE : R_LUCK);
        first = false;

        moveSelectColl:

        Move move = config.manualMoves || config.logMode == LOG_READ ? chooseMove(config) :
                  (config.competitive ? getCompetitiveCollMove(states, num, collected, config) : bestMove);

        if (!config.manualMoves && config.verbose) std::cout << "Move: " << move.toString() << std::endl;

        switch (move.id)
        {
        case M_LIST_OPTIONS:
        {
            std::unordered_set<int> mentioned;
            for (Move& option : getCollMoveOptions(state.goods, left))
            {
                if (mentioned.insert(option.id).second) std::cout << option.toString(true) << std::endl;
            }
            fail = true;
            break;
        }

        case M_STOP_COLL:
            if (config.logMode == LOG_WRITE) config.logOut << move.toString() << std::endl;

            if (config.evalLM) state.updateExpScore(collected * (num + 1) + getScore(state.free, state.goods), R_MISTAKE);
            break;
        
        case M_CONT_COLL:
            if (state.goods > 0 && left > 0)
            {
                if (config.logMode == LOG_WRITE) config.logOut << move.toString() << std::endl;

                if (config.evalLM) state.updateExpScore(collected * (num + 1) + getCollContScore(state.free, state.goods, num, left), R_MISTAKE);
                cont = true;
                --state.goods;
                int newHits;
                if (config.manualRolls || config.logMode == LOG_READ) newHits = chooseNumNewHits(num, left, config);
                else newHits = randNumNewHits(left);

                if (config.logMode == LOG_WRITE) config.logOut << newHits << std::endl;

                if (!config.manualRolls && config.verbose)
                {
                    std::cout << "Number of " << numNames[num] << " rolled: " << newHits << std::endl;
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
            if (config.verbose && move.id != M_LIST_OPTIONS) std::cout << "Invalid move." << std::endl;
            fail = false;
            goto moveSelectColl;
        }
    }

    assert(!fail);

    int points = collected * (num + 1);
    state.score += points;

    if (config.verbose)
    {
        std::cout << "Won " << points << " points." << std::endl;
    }
}

void simMove(std::vector<State>& states, int code, int extraDice, int points, Config& config)
{
    State& state = states[config.player];

    bool instaDone = isDone(code);
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
        rolls = randNumRolls(code, extraDice, state.goods);
        won = isDone(code);
    }

    if (!instaDone && state.goods > 0 && config.logMode == LOG_WRITE) config.logOut << (won ? rolls : -1) << std::endl;

    if (!instaDone && !config.manualRolls && config.verbose)
    {
        if (won) std::cout << "Took " << rolls << " rolls to complete the combination. " << std::endl;
        else std::cout << "Did not complete the combination." << std::endl;
    }

    state.goods -= rolls;
    if (won) state.score += points;

    if (config.evalLM) state.updateExpScore(getScore(state.free, state.goods), instaDone ? R_NONE : R_LUCK);

    if (won && config.verbose)
    {
        std::cout << "Won " << points << " points." << std::endl;
    }
}

void simTurn(std::vector<State>& states, Config& config)
{
    State& state = states[config.player];
    int diceOrdCode;

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
        std::cout << "Expected final score: " << state.expScore << std::endl;
    }

    if (!config.manualMoves && config.competitive) findOthersCumDistr(states, config);

    state.goods += GOODS_PER_TURN;

    diceRoll:

    if (config.manualRolls || config.logMode == LOG_READ) diceOrdCode = chooseOrdCode(config);
    else diceOrdCode = randOrdCode();

    if (config.logMode == LOG_WRITE) printOrdCode(diceOrdCode, config, true);
    if (!config.manualRolls && config.verbose) printOrdCode(diceOrdCode, config, false);

    Move bestMove = getMove(state.free, state.goods, diceOrdCode);

    if (config.evalLM) state.updateExpScore(bestMove.score, R_LUCK);

    if (config.verbose)
    {
        std::cout << "Goods: " << state.goods << std::endl;
    }

    bool fail = false;

    moveSelect:

    Move move = config.manualMoves || config.logMode == LOG_READ ? chooseMove(config) :
        (config.competitive ? getCompetitiveMove(states, diceOrdCode, config) : bestMove);

    if (!config.manualMoves && config.verbose) std::cout << "Move: " << move.toString() << std::endl;

    switch (move.id)
    {
    case M_LIST_OPTIONS:
    {
        std::unordered_set<int> mentioned;
        for (Move& option : getMoveOptions(state.free, state.goods, diceOrdCode))
        {
            if (mentioned.insert(option.id).second) std::cout << option.toString(true) << std::endl;
        }
        fail = true;
        break;
    }

    case M_ERROR:
    case M_STOP_COLL:
    case M_CONT_COLL:
        fail = true;
        break;

    case M_REROLL:
        if (state.goods > 0 && !MISERE)
        {
            if (config.logMode == LOG_WRITE) config.logOut << move.toString() << std::endl;

            --state.goods;
            if (config.evalLM) state.updateExpScore(getContScore(state.free, state.goods), R_MISTAKE, bestMove.toString());
            goto diceRoll;
        }
        else fail = true;
        break;
    
    default:
        if (move.id >= 0 && move.id < NUM_COMBOS && isFree(state.free, move.id))
        {
            if (config.logMode == LOG_WRITE) config.logOut << move.toString() << std::endl;

            state.free = setUsed(state.free, move.id);
            int num = combos[move.id]->getCollectNumber();
            if (num != NONUM)
            {
                    int collected = ordCodeOccurs[diceOrdCode][num];
                if (config.evalLM) state.updateExpScore(collected * (num + 1) + getCollScore(state.free, state.goods, num, NUM_DICE - collected), R_MISTAKE, bestMove.toString());
                simColl(states, num, collected, config);
            }
            else
            {
                int newCode = ordCodeRemOrdCode[move.target.ordCode][diceOrdCode];
                int extraDice = NUM_DICE - diceInOrdCode[move.target.ordCode];
                if (config.evalLM) state.updateExpScore(getMoveScore(state.free, state.goods, newCode, extraDice, move.target.points), R_MISTAKE, bestMove.toString());
                simMove(states, newCode, extraDice, move.target.points, config);
            }
        }
        else fail = true;
    }

    if (config.manualMoves && fail)
    {
        if (config.verbose && move.id == M_LIST_OPTIONS) std::cout << "Invalid move" << std::endl;
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

    file.setf(std::ios::scientific);
    file.precision(FILE_PRECISION);

    std::cout << "Storing model in " << name << "." << std::endl;

    for (int i = 0; i <= MAX_EXTRA_DICE; ++i)
    {
        for (int j = 0; j < VALID_CODES; ++j)
        {
            if (i + diceInCode[j] > NUM_DICE) continue;
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

    std::cout << "Stored model." << std::endl;
}

bool loadModel()
{
    std::string name = modelName();
    std::ifstream file(name.c_str());

    if (!file) return false;

    std::cout << "Loading model from " << name << "." << std::endl;

    for (int i = 0; i <= MAX_EXTRA_DICE; ++i)
    {
        for (int j = 0; j < VALID_CODES; ++j)
        {
            if (i + diceInCode[j] > NUM_DICE) continue;
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
            std::cout << "Mean: " << stats.mean << std::endl;
            std::cout << "Stdev: " << stats.stdev << std::endl;

            std::cout.precision(useRank ? 1 : 0);
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
            std::cout.precision(IO_PRECISION);

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
            std::cout << "Expected score: " << getInitialScore() << std::endl;
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

    std::cout.setf(std::ios::fixed);
    std::cout.precision(IO_PRECISION);

    std::cout << "Normal play or misere play (0 / 1): ";
    std::cin >> MISERE;

    std::cout << "Number of trails per state or 0 for exact model (which costs " << VALID_ORD_CODES - START_FULL_ORD_CODES << " trials per state): ";
    std::cin >> NUM_TRIALS;

    genCodes();
    setCombos();

    if (!loadModel())
    {
        computeModel();
        storeModel();
    }

    shell();

    return 0;
}
