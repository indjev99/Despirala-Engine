#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <math.h>
#include <stdlib.h>
#include <time.h>

const int NUM_DICE = 6;
const int NUM_SIDES = 6;
const int MAX_ROLL = 46656; // NUM_DICE ^ NUM_SIDES
const int MAX_CODE = 1 << (NUM_SIDES + NUM_DICE);
const int DIFF_CODES = 30;
const int NUM_COMBOS = 14;
const int NUM_MASKS = 1 << NUM_COMBOS;
const int GOODS_PER_TURN = 5;
const int POINTS_PER_GOOD = 1;
const int MAX_GOODS = NUM_COMBOS * GOODS_PER_TURN;
const int NUM_TRIALS_1 = 100000;
int NUM_TRIALS_2;

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
    virtual int getNumber() const
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
    int getNumber() const
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
        for (int i = 0; i < args.size(); ++i)
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

void randOcc(Occurs& occurs)
{
    occurs.fill(0);
    for (int i = 0; i < NUM_DICE; ++i)
    {
        ++occurs[rand() % NUM_SIDES];
    }
}

bool randRemOcc(int extraDice, Occurs& occurs)
{
    int numFreeDice = std::accumulate(occurs.begin(), occurs.end(), extraDice);
    bool changed = false;
    for (int i = 0; i < numFreeDice; ++i)
    {
        int side = rand() % NUM_SIDES;
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

const int MAX_EXTRA_DICE = 2;
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

#define STOP_COLL -1
#define CONT_COLL -2
#define REROLL -3
#define ERROR -100

struct Move
{
    int id;
    std::vector<int> args;
    double score;

    Move(int id, double score):
        id(id),
        args({}),
        score(score) {}
    Move(int id, const std::vector<int>& args, double score):
        id(id),
        args(args),
        score(score) {}
    std::string toString() const
    {
        std::string s;
        switch (id)
        {
        case ERROR:
            s = "Error";
            break;
        case STOP_COLL:
            s = "Stop collecting";
            break;
        case CONT_COLL:
            s = "Continue collecting";
            break;
        case REROLL:
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
};

bool operator<(const Move& a, const Move& b)
{
    return a.score < b.score;
}

double getScore(int free, int goods);
double getScoreCont(int free, int goods);
double getScoreColl(int free, int goods, int num, int left);

Move getMoveColl(int free, int goods, int num, int left)
{
    Move best = Move(STOP_COLL, getScore(free, goods));
    if (goods && left)
    {
        double sc = 0;
        for (int i = 0 ; i <= left; ++i)
        {
            sc += leftDistr[left][i] * ((left - i) * num + getScoreColl(free, goods - 1, num, i));
        }
        Move option = Move(CONT_COLL, sc);
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

double getMoveScore(int free, int goods, const Occurs& occurs, int extraDice, double reward)
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
    Move best = Move(ERROR, 0);
    for (int i = 0; i < NUM_COMBOS; ++i)
    {
        if (!isFree(free, i)) continue;
        int newFree = setUsed(free, i);
        int num = combos[i]->getNumber();
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
        Move option = Move(REROLL, getScoreCont(free, goods - 1));
        best = std::max(best, option);
    }
    return best;
}

int cnt = 0;

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
    ++cnt;
    //std::cerr << "score[" << free << "][" << goods << "] = " << score[free][goods] << "; // " << cnt << "\n";
    //if (cnt % 10000 == 0) std::cerr << cnt << "\n";
    return score[free][goods];
}

double getScore(int free, int goods)
{
    return free ? getScoreCont(free, goods + GOODS_PER_TURN) : getScoreCont(free, goods);
}

void findExpectedScores()
{
    getScore(NUM_MASKS - 1, 0);
}

struct State
{
    int free;
    int goods;
    int score;
    bool fail;

    State():
        free(NUM_MASKS - 1),
        goods(0),
        score(0),
        fail(false) {}
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

const std::string numNames[NUM_SIDES] = {"ones", "twos", "threes", "fours", "fives", "sixes"};
void simCollMoves(State& s, int num, int left, int verbosity)
{
    bool cont = true;
    while (cont)
    {
        cont = false;
        Move mv = getMoveColl(s.free, s.goods, num, left);

        if (verbosity >= 4)
        {
            std::cout << "Have " << NUM_DICE - left << " " << numNames[num - 1] << " and " << s.goods << " goods." << std::endl;
        }
        if (verbosity >= 3)
        {
            std::cout << "Move: " << mv.toString() << std::endl;
        }

        switch (mv.id)
        {
        case STOP_COLL:
            break;
        case CONT_COLL:
            if (s.goods)
            {
                cont = true;
                --s.goods;
                int newHits = 0;
                for (int i = 0; i < left; ++i)
                {
                    if (rand() % NUM_SIDES == 0) ++newHits;
                }
                left -= newHits;
            }
            else s.fail = true;
            break;
        default:
            s.fail = true;
        }
    }
    int reward = num * (NUM_DICE - left);
    s.score += reward;

    if (verbosity >= 3)
    {
        std::cout << "Won " << reward << " points." << std::endl;
    }
}

bool isDone(const Occurs& occurs)
{
    return std::all_of(occurs.begin(), occurs.end(), [](int n) {return n == 0;});
}

void simRegMove(State& s, Occurs& occurs, int ed, int reward, int verbosity)
{
    int rolls = 0;
    while (!isDone(occurs) && rolls < s.goods)
    {
        if (verbosity >= 5)
        {
            std::cout << "Need: ";
            printOccurs(occurs);
        }

        ++rolls;
        randRemOcc(ed, occurs);
    }
    s.goods -= rolls;
    if (isDone(occurs)) s.score += reward;

    if (verbosity >= 3)
    {
        std::cout << "Took " << rolls << " rolls. " << std::endl;
        if (isDone(occurs)) std::cout << "Won " << reward << " points." << std::endl;
        else std::cout << "Didn't complete the combination." << std::endl;
    }
}

int simGame(int verbosity=-1)
{
    State s;
    Occurs diceOccurs;
    bool finishedTurn = true;

    while (!s.fail && s.free)
    {
        if (finishedTurn) s.goods += GOODS_PER_TURN;
        randOcc(diceOccurs);
        Move mv = getMove(s.free, s.goods, diceOccurs);

        if (verbosity >= 1)
        {
            std::cout << "Have " << s.score << " points and " << s.goods << " goods." << std::endl;
            if (verbosity >= 2)
            {
                std::cout << "Rolled: ";
                printOccurs(diceOccurs);
            }
            std::cout << "Move: " << mv.toString() << std::endl;
        }

        int id = mv.id;
        switch (id)
        {
        case ERROR:
        case STOP_COLL:
        case CONT_COLL:
            s.fail = true;
            break;
        case REROLL:
            if (s.goods)
            {
                --s.goods;
                finishedTurn = false;
            }
            else s.fail = true;
            break;
        default:
            if (isFree(s.free, id))
            {
                s.free = setUsed(s.free, id);
                int num = combos[id]->getNumber();
                if (num)
                {
                    simCollMoves(s, num, NUM_DICE - diceOccurs[num - 1], verbosity);
                }
                else
                {
                    Target t = combos[id]->getTargetByArgs(mv.args);
                    Occurs newOccurs = t.occurs;
                    int extraDice = calcExtraDice(newOccurs);
                    remOcc(diceOccurs, newOccurs);
                    simRegMove(s, newOccurs, extraDice, t.points, verbosity);
                }
                finishedTurn = true;
            }
            else s.fail = true;

            if (verbosity >= 1 && finishedTurn) std::cout << std::endl;
        }
    }

    int score = s.score + s.goods * POINTS_PER_GOOD;

    if (verbosity >= 0) std::cout << "Final score: " << score << std::endl;

    return score;
}

#define INF 1e9

struct Stats
{
    double mean;
    double stdev;
    int median;
    int perc5 = INF;
    int perc95 = -INF;
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

    std::cout << "Loading the model from " << name << "." << std::endl;

    if (!file) return false;

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

void trainModel()
{
    std::cout << "Training the model." << std::endl;

    generateLefts();
    std::cout << "Generated possible things left to roll." << std::endl;

    findRollsDistr();
    std::cout << "Found the distribution of dice rolls needed for them." << std::endl;

    findLeftDistr();
    std::cout << "Found the distribution of number of dice left after roll." << std::endl;

    findExpectedScores();
    std::cout << "Found the expected scores." << std::endl;
    std::cout << "Considered " << cnt << " cases." << std::endl;
}

const std::string help = "Possible commands: stop, help, example, test, expected.";
void shell()
{
    std::string cmd;
    std::cout << "\n" << help << std::endl;
    while (true)
    {
        std::cout << "\nEnter command: ";
        std::cin >> cmd;

        if (cmd == "stop") break;
        if (cmd == "help")
        {
            std::cout << help << std::endl;
        }
        else if (cmd == "example")
        {
            int verbosity;
            std::cout << "Verbosity (0 - 5): ";
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
    }
}

int main()
{
    srand(time(0));

    std::cout << "Number of trails per state, higher leads to a better model: ";
    std::cin >> NUM_TRIALS_2;

    if (!loadModel())
    {
        trainModel();
        storeModel();
    }

    shell();

    return 0;
}
