#include <iostream>
#include <fstream>
#include <sstream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <cstdio>
#include <memory>
#include <random>
#include <list>
#include <string>
#include <algorithm>
#include <iomanip>
#include <chrono>

using namespace std;

typedef long long int i64;

vector<i64> BestSuite;

using namespace std;
struct graph_data {
    constexpr static i64 MaxDegree = 5;
    i64 V, E;
    std::map<i64, std::map<i64, i64>> edges;
    graph_data(std::istream& src) {
        src >> V >> E;
        for (i64 i = 0; i < E; ++i) {
            i64 u, v, d;
            src >> u >> v >> d;
            --u, --v;
            edges[u][v] = d;
            edges[v][u] = d;
        }
        for ([[maybe_unused]] auto& [u, connects] : edges) {
        }
    }

    void dump(ostream& ss) {
        ss << V << " " << E << endl;
        map<pair<i64, i64>, i64> q;
        for (auto& [u, t] : edges) {
            for (auto& [v, d] : t) {
                q[{min(u, v), max(u, v)}] = d;
            }
        }
        for (auto [t, d] : q) {
            auto [u, v] = t;
            ss << u << " " << v << " " << d << endl;
        }
    }
};

struct grid_data {
    i64 DayType = 0;
    i64 N_div = 0;
    i64 N_pattern = 0;
    i64 sigma_ele = 0;
    i64 Delta_event = 0;
    double p_event = 0.0;
    std::vector<std::vector<i64>> pw_predict;
    i64 N_grid, C_grid_init, C_grid_max, V_grid_max;
    std::vector<i64> x, pattern;
    grid_data() = default;
    grid_data(const grid_data&) = default;
    grid_data(std::istream& src) {
        src >> DayType;
        src >> N_div >> N_pattern >> sigma_ele >> p_event >> Delta_event;
        pw_predict.resize(N_pattern);
        for (i64 i = 0; i < N_pattern; ++i) {
            pw_predict[i].resize(N_div);
            for (i64 j = 0; j < N_div; ++j) {
                src >> pw_predict[i][j];
            }
        }
        src >> N_grid >> C_grid_init >> C_grid_max >> V_grid_max;
        x.resize(N_grid);
        pattern.resize(N_grid);
        for (i64 i = 0; i < N_grid; ++i) {
            src >> x[i] >> pattern[i];
            --x[i];
            --pattern[i];
        }
    }
    void dump(ostream& ss) {
        ss << DayType << endl;
        ss << N_div << " " << N_pattern << " " << sigma_ele << " " << p_event << " " << Delta_event << endl;
        for (i64 i = 0; i < N_pattern; ++i) {
            for (i64 j = 0; j < N_div; ++j) {
                ss << pw_predict[i][j] << " ";
            }
            ss << endl;
        }
        ss << N_grid << " " << C_grid_init << " " << C_grid_max << " " << V_grid_max << endl;
        for (i64 i = 0; i < N_grid; ++i) {
            ss << x[i] << " " << pattern[i] << endl;
        }

    }
    grid_data& operator=(const grid_data&) = default;
};

struct EV_data {
    i64 N_EV = 0;
    i64 C_EV_init = 0;
    i64 C_EV_max = 0;
    i64 V_EV_max = 0;
    i64 N_trans_max = 0;
    i64 Delta_EV_move = 0;
    std::vector<i64> pos;
    EV_data() = default;
    EV_data(std::istream& src, bool isA) {
        if (isA) {
            src >> N_EV >> C_EV_init >> C_EV_max >> V_EV_max >> Delta_EV_move;
        }
        else {
            src >> N_EV >> C_EV_init >> C_EV_max >> V_EV_max >> N_trans_max >> Delta_EV_move;
        }
        pos.resize(N_EV);
        for (i64 i = 0; i < N_EV; ++i) {
            src >> pos[i];
            --pos[i];
        }
    }

    void dump(ostream& ss, bool isA) {
        if (isA) {
            ss << N_EV << " " << C_EV_init << " " << C_EV_max << " " << V_EV_max << " " << Delta_EV_move << endl;
        }
        else {
            ss << N_EV << " " << C_EV_init << " " << C_EV_max << " " << V_EV_max << " " << N_trans_max << " " << Delta_EV_move << endl;
        }
        for (i64 i = 0; i < N_EV; ++i) {
            ss << pos[i] << endl;
        }
    }
};

struct B {
    graph_data graph;
    grid_data grid;
    EV_data EV;

    double p_const_trans = 0.0;
    i64 T_last = 0;
    i64 P_trans = 0;
    double gamma = 0.0;
    i64 S_ref_ele = 0;
    i64 S_ref_trans = 0;
    i64 T_max = 0;

    B(std::istream& src, bool isA = false) :graph(src), grid(src), EV(src, isA) {
        if (isA) {
            src >> gamma;
            src >> T_max;
        }
        else {
            src >> p_const_trans >> T_last;
            src >> P_trans >> gamma >> S_ref_ele >> S_ref_trans;
            src >> T_max;
        }
    }

    void dump(ostream& ss, bool isA) {
        graph.dump(ss);
        grid.dump(ss);
        EV.dump(ss, isA);

        if (isA) {
            ss << gamma << endl;
            ss << T_max << endl;
        }
        else {
            ss << p_const_trans << " " << T_last << endl;
            ss << P_trans << " " << gamma << " " << S_ref_ele << " " << S_ref_trans << endl;
            ss << T_max << endl;
        }
    }
};

struct carinfo {
    i64 charge = 0;
    i64 u = 0;
    i64 v = 0;
    i64 dist_from_u = 0;
    i64 dist_to_v = 0;
    i64 N_adj = 0;
    std::vector<i64> a;
    i64 N_order = 0;
    std::vector<i64> o;

    void load(std::istream& src, bool isA, [[maybe_unused]] i64 C_EV_max = 25000, [[maybe_unused]] i64 V = 225, [[maybe_unused]] i64 MaxDegree = 5, [[maybe_unused]] i64 N_trans_max = 4, [[maybe_unused]] i64 T_last = 900) {
        src >> charge;
        src >> u >> v >> dist_from_u >> dist_to_v;
        --u, --v;
        src >> N_adj; a.resize(N_adj);
        for (i64 i = 0; i < N_adj; ++i) {
            src >> a[i];
            --a[i];
        }
        if (!isA) {
            src >> N_order; o.resize(N_order);
            for (i64 i = 0; i < N_order; ++i) {
                src >> o[i];
                --o[i];
            }
        }
    }

    void dump(ostream& ss, bool isA) {
        ss << charge << endl;
        ss << u << " " << v << " " << dist_from_u << " " << dist_to_v << endl;

        ss << N_adj << endl;
        for (i64 i = 0; i < N_adj; ++i) {
            ss << a[i] << " ";
        }
        ss << endl;

        if (!isA) {
            ss << N_order << endl;
            for (i64 i = 0; i < N_order; ++i) {
                ss << o[i] << " ";
            }
            ss << endl;
        }
    }
};

struct grid_info {
    i64 N_grid = 0;
    std::vector<i64> x, y;
    std::vector<i64> pw_actual;
    std::vector<i64> pw_excess, pw_buy;
    grid_info() = default;
    grid_info(i64 N_grid)
        :N_grid(N_grid), x(N_grid), y(N_grid), pw_actual(N_grid), pw_excess(N_grid), pw_buy(N_grid) {}
    void load(std::istream& src, [[maybe_unused]] i64 V = 225, [[maybe_unused]] i64 C_grid_max = 50000) {
        for (i64 i = 0; i < N_grid; ++i) {
            src >> x[i] >> y[i] >> pw_actual[i] >> pw_excess[i] >> pw_buy[i];
            --x[i];
        }
    }

    void dump(ostream& ss) {
        for (i64 i = 0; i < N_grid; ++i) {
            ss << x[i] << " " << y[i] << " " << pw_actual[i] << " " << pw_excess[i] << " " << pw_buy[i] << endl;
        }
    }
};

struct EV_info {
    i64 N_EV = 0;
    std::vector<carinfo> c;
    EV_info() = default;
    EV_info(i64 N_EV)
        :N_EV(N_EV), c(N_EV) {}
    void load(std::istream& src, bool isA) {
        for (i64 i = 0; i < N_EV; ++i) {
            c[i].load(src, isA);
        }
    }
    void dump(ostream& ss, bool isA) {
        for (i64 i = 0; i < N_EV; ++i) {
            c[i].dump(ss, isA);
        }
    }
};

struct order_info {
    i64 N_order = 0;
    std::vector<i64> id;
    std::vector<i64> w;
    std::vector<i64> z;
    std::vector<i64> state;
    std::vector<i64> time;

    void load(std::istream& src, [[maybe_unused]] i64 V = 225, [[maybe_unused]] i64 T_last = 900) {
        src >> N_order;
        id.resize(N_order);
        w.resize(N_order);
        z.resize(N_order);
        state.resize(N_order);
        time.resize(N_order);
        for (i64 i = 0; i < N_order; ++i) {
            src >> id[i] >> w[i] >> z[i] >> state[i] >> time[i];
            --w[i], --z[i];
        }
    }

    void dump(ostream& ss) {
        ss << N_order << endl;
        for (i64 i = 0; i < N_order; ++i) {
            ss << id[i] << " " << w[i] << " " << z[i] << " " << state[i] << " " << time[i] << endl;
        }
    }
};

struct graph_summary {
    vector<vector<i64>> len;
    vector<vector<i64>> next;
    vector<i64> nanogrid_pos;
    i64 diameter = 0;
    i64 cover_radius = 0;
    graph_summary(const graph_data& graph, const grid_data& grid) :
        len(graph.V, std::vector<i64>(graph.V, (i64)1e9)),
        next(graph.V, std::vector<i64>(graph.V)),
        nanogrid_pos(grid.N_grid) {
        const i64 V = graph.V;
        for (i64 i = 0; i < V; ++i)
            len[i][i] = 0;
        for (i64 i = 0; i < V; ++i)
            for (i64 j = 0; j < V; ++j)
                next[i][j] = j;
        for (const auto& [u, u_edges] : graph.edges)
            for (const auto& [v, length] : u_edges) {
                len[u][v] = length;
                len[v][u] = length;
            }
        for (i64 k = 0; k < V; ++k)
            for (i64 i = 0; i < V; ++i)
                for (i64 j = 0; j < V; ++j)
                    if (len[i][j] > len[i][k] + len[k][j]) {
                        len[i][j] = len[i][k] + len[k][j];
                        next[i][j] = next[i][k];
                    }
        nanogrid_pos = grid.x;
        for (i64 i = 0; i < V; ++i)
            for (i64 j = 0; j < V; ++j)
                diameter = max(len[i][j], diameter);
        for (i64 i = 0; i < V; ++i) {
            i64 min_len = (i64)1e9;
            for (i64 j = 0; j < (i64)nanogrid_pos.size(); ++j)
                min_len = min(min_len, len[i][j]);
            cover_radius = max(cover_radius, min_len);
        }
    }
};

struct strategy : public B {
    const graph_summary& gs;
    vector<list<string>> command_queue;
    strategy(const B& p, const graph_summary& gs) : B(p), gs(gs),
        command_queue(B::EV.N_EV) {}
    virtual void command(const grid_info& g_i, const EV_info& ev_i, const order_info& order_i) = 0;
    virtual void initialize() {
        for (auto& queue : command_queue) queue.clear();
    }
    bool is_free(i64 EV_index) {
        if (command_queue[EV_index].size() > 0) {
            return false;
        }
        return true;
    }
    string dequeue(const EV_info& ev_i) {
        string ret = "";
        for (i64 i = 0; i < ev_i.N_EV; ++i)
            ret += dequeue(i) + "\n";
        return ret;
    }
    string dequeue(i64 EV_index) {
        string ret;
        if (command_queue[EV_index].size() > 0) {
            ret = command_queue[EV_index].front();
            command_queue[EV_index].pop_front();
        }
        else {
            ret = "stay";
        }
        return ret;
    }
    void enqueue(i64 EV_index, const string& cmd) {
        command_queue[EV_index].push_back(cmd);
    }
    void enqueue(i64 EV_index, const string& cmd, i64 repeat) {
        for (i64 i = 0; i < repeat; ++i)
            command_queue[EV_index].push_back(cmd);
    }
    void enqueue(i64 EV_index, list<string>&& cmd_list) {
        command_queue[EV_index].splice(command_queue[EV_index].end(), cmd_list);
    }
};

vector<string> split_command(const string& command_pack) {
    vector<string> ret;
    stringstream reader(command_pack);
    string line;
    while (getline(reader, line)) {
        if (line == "") continue;
        else ret.emplace_back(line);
    }
    return ret;
}

struct TStrategy : public strategy {
public:
    bool isA = false;

public:
    map<i64, i64> Orders;
    map<i64, i64> Grids;

public:
    set<pair<i64, i64>> EvsQueue;
    set<pair<i64, i64>> EvsChargeQueue;

public:
    set<i64> AssignedOrders;

public:
    vector<set<i64>> EvTargetOrder;
    vector<set<i64>> EvCarringOrder;

public:
    vector<i64> EvTargetGrid;
    vector<i64> EvChargeGrid;

public:
    i64 tt = 0;
    i64 tm = 0;

public:
    i64 EvLimit = 0;

public:
    vector<vector<i64>> ConsumptionByGrid;
    vector<map<i64, map<i64, i64>>> ConsumptionByEv;

public:
    vector<i64> GridSumBuy;

public:
    vector<vector<pair<pair<i64, i64>, i64>>> Incidents;

public:
    vector<pair<i64, i64>> GridLimits;
    vector<i64> GridNextLoss;
    vector<i64> GridChargeNeeded;

public:
    vector<pair<double, double>> Scores;
    pair<double, double> PredictedScore;

public:
    set<i64> WeatherChanges;

public:
    TStrategy(const B& b, const graph_summary& gs, i64 qt, bool _isA, vector<pair<double, double>>& scores, pair<double, double> predictedScore) : strategy(b, gs), tt(qt), isA(_isA), Scores(scores), PredictedScore(predictedScore) {
        ConsumptionByGrid = vector<vector<i64>>(b.grid.N_grid, vector<i64>(T_max));
        ConsumptionByEv = vector<map<i64, map<i64, i64>>>(EV.N_EV);
    }

    void initialize() {
        strategy::initialize();

        EvLimit = GetEvLimit();

        AssignedOrders.clear();
        EvTargetOrder = vector<set<i64>>(EV.N_EV);
        EvCarringOrder = vector<set<i64>>(EV.N_EV);
        EvTargetGrid = vector<i64>(EV.N_EV, -1);
        EvChargeGrid = vector<i64>(EV.N_EV, -1);

        GridSumBuy = vector<i64>(grid.N_grid);

        Incidents = vector<vector<pair<pair<i64, i64>, i64>>>(grid.N_grid);

        GridLimits = vector<pair<i64, i64>>(grid.N_grid);
        GridNextLoss = vector<i64>(grid.N_grid);
        GridChargeNeeded = vector<i64>(grid.N_grid);
    }

    void SimulateCharging(i64 gridId, i64 startTm, i64 chargeNeeded, bool rollback, i64 evId = -1) {
        if (chargeNeeded < 0) {
            chargeNeeded *= -1;
            rollback = !rollback;
        }

        for (i64 ftm = startTm; ftm < T_max; ftm++) {
            i64 charge = min(chargeNeeded, EV.V_EV_max);
            if (charge <= 0) {
                break;
            }
            chargeNeeded -= charge;

            if (evId >= 0) {
                if (!rollback) {
                    ConsumptionByEv[evId][gridId][ftm] += charge;
                }
                else {
                    ConsumptionByEv[evId][gridId][ftm] -= charge;
                }
                if (!ConsumptionByEv[evId][gridId][ftm]) {
                    ConsumptionByEv[evId][gridId].erase(ftm);
                }
            }

            ConsumptionByGrid[gridId][ftm] += (rollback ? -1 : 1) * charge;
        }
    }

    i64 GetPw(const grid_info& grid_i, i64 gridId, i64 tm, i64 curTm) {
        i64 patternId = grid.pattern[gridId];
        i64 div = T_max / grid.N_div;

        if (curTm && tm / div == (curTm - 1) / div && abs(grid.pw_predict[patternId][tm / div] - grid_i.pw_actual[gridId]) > 500) {
            return grid_i.pw_actual[gridId];
        }
        return grid.pw_predict[patternId][tm / div];
    }


    pair<i64, i64> CalculateLosses(const grid_info& grid_i, i64 gridId, i64 curTm, i64 deltaTm = 0, i64* charge = nullptr, i64* minCharge = nullptr, i64* maxCharge = nullptr, bool allLosses = false) {
        i64 sumExcess = 0;
        i64 sumBuy = 0;

        i64 c = grid_i.y[gridId];

        if (minCharge) {
            *minCharge = 10 * grid.C_grid_max;
        }

        if (maxCharge) {
            *maxCharge = 0;
        }


        for (i64 tm = curTm; tm < T_max; tm++) {
            i64 delta = GetPw(grid_i, gridId, tm, curTm) - ConsumptionByGrid[gridId][tm];

            i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - c);
            if (excess < 0) {
                excess = 0;
            }

            i64 buy = -delta - min(grid.V_grid_max, c);
            if (buy < 0) {
                buy = 0;
            }

            c += delta;
            if (c < 0) {
                c = 0;
            }
            if (c > grid.C_grid_max) {
                c = grid.C_grid_max;
            }

            if (tm - curTm >= deltaTm) {
                sumExcess += excess;
                sumBuy += buy;

                if (minCharge) {
                    *minCharge = min(*minCharge, c);
                }

                if (maxCharge) {
                    *maxCharge = max(*maxCharge, c);
                }
            }
            else if (allLosses) {
                sumExcess += excess;
                sumBuy += buy;
            }
        }

        if (charge) {
            *charge = c;
        }

        return { sumExcess, sumBuy };
    }

    i64 SimulateDelivery(const EV_info& ev_i, const order_info& order_i, map<i64, i64>& orders, i64 evId) {
        auto& ev = ev_i.c[evId];

        i64 shift = ev_i.c[evId].u == ev_i.c[evId].v ? 0 : (ev_i.c[evId].u == EvTargetGrid[evId] ? ev_i.c[evId].dist_from_u : ev_i.c[evId].dist_to_v);
        i64 x = ev_i.c[evId].u == ev_i.c[evId].v ? ev_i.c[evId].u : EvTargetGrid[evId];

        i64 res = shift;

        auto targetOrder = EvTargetOrder[evId];
        auto carringOrder = EvCarringOrder[evId];

        while (!targetOrder.empty() || !carringOrder.empty()) {
            map<i64, vector<i64>> dest;
            for (i64 orderId : targetOrder) {
                dest[order_i.w[orders[orderId]]].push_back(orderId);
            }
            for (i64 orderId : carringOrder) {
                dest[order_i.z[orders[orderId]]].push_back(orderId);
            }


            i64 closestDist = (i64)1e18;
            i64 closestGrid = -1;

            for (auto [grid, _] : dest) {
                i64 dist = gs.len[x][grid];
                if (closestDist > dist) {
                    closestDist = dist;
                    closestGrid = grid;
                }
            }

            if (closestGrid < 0) {
                cerr << "closestGrid < 0" << endl;
                throw 1;
            }
            res += gs.len[x][closestGrid];
            x = closestGrid;

            for (auto orderId : dest[closestGrid]) {
                if (order_i.w[orders[orderId]] == closestGrid) {
                    targetOrder.erase(orderId);
                    carringOrder.insert(orderId);
                    res += 1;
                    continue;
                }
                if (order_i.z[orders[orderId]] == closestGrid) {
                    carringOrder.erase(orderId);
                    continue;
                }

                cerr << "SimulateDelivery" << endl;
                throw 1;
            }
        }

        return res;
    }

    i64 EsimateDeliveryCost(const EV_info& ev_i, const order_info& order_i, map<i64, i64>& orders, i64 evId, i64 orderId) {
        i64 now = SimulateDelivery(ev_i, order_i, orders, evId);
        EvTargetOrder[evId].insert(orderId);
        i64 nxt = SimulateDelivery(ev_i, order_i, orders, evId);
        EvTargetOrder[evId].erase(orderId);

        if (nxt >= T_max - tm) {
            return 1'000'000;
        }

        return nxt - now;
    }

private:
    void PreUpdate(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            GridSumBuy[gridId] += grid_i.pw_buy[gridId];
        }
    }

private:
    i64 GetEvLimit() {
        if (isA) {
            return 0;
        }

        if (tt < 0) {
            return -tt;
        }

        return BestSuite[tt];
    }

private:
    map<i64, i64> GetOrdersMapping(const order_info& order_i) {
        map<i64, i64> orders;
        for (i64 i = 0; i < order_i.N_order; i++) {
            orders[order_i.id[i]] = i;
        }
        return orders;
    }

    map<i64, i64> GetGridsMapping(const grid_info& grid_i) {
        map<i64, i64> grids;
        for (i64 i = 0; i < grid_i.N_grid; i++) {
            grids[grid_i.x[i]] = i;
        }

        return grids;
    }

private:
    void ConstructEvsQueues(const EV_info& ev_i) {
        if (!tm) {
            for (i64 evId = 0; evId < EV.N_EV; evId++) {
                if (EvChargeGrid[evId] < 0) {
                    continue;
                }

                EvsQueue.insert({ ev_i.c[evId].charge, evId });
            }

            for (i64 evId = 0; evId < EV.N_EV && (i64)EvsQueue.size() < GetEvLimit(); evId++) {
                if (EvChargeGrid[evId] >= 0) {
                    continue;
                }
                EvsQueue.insert({ ev_i.c[evId].charge, evId });
            }

            for (i64 evId = 0; evId < EV.N_EV; evId++) {
                if (EvsQueue.count({ ev_i.c[evId].charge, evId })) {
                    continue;
                }
                EvsChargeQueue.insert({ ev_i.c[evId].charge, evId });
            }

            return;
        }

        if (grid.DayType == 1 && !EvsQueue.empty() && PredictedScore != make_pair(0.0, 0.0)) {
            i64 sumBuy = 0;
            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                sumBuy += GridSumBuy[gridId];
            }

            bool ok = true;
            for (auto [trans, ele] : Scores) {
                if (trans > PredictedScore.first && ele > PredictedScore.second - sumBuy - 3 * 50 * 1'000 * WeatherChanges.size()) {
                    ok = false;
                    break;
                }
            }

            if (!ok) {
                cerr << tm << ": " << "!ok" << endl;

                for (auto t : EvsQueue) {
                    EvsChargeQueue.insert(t);
                }

                EvsQueue.clear();
            }
        }

        set<pair<i64, i64>> evsQueue;
        for (auto [_, evId] : EvsQueue) {
            evsQueue.insert({ ev_i.c[evId].charge, evId });
        }

        set<pair<i64, i64>> evsChargeQueue;
        for (auto [_, evId] : EvsChargeQueue) {
            evsChargeQueue.insert({ ev_i.c[evId].charge, evId });
        }

        EvsQueue = evsQueue;
        EvsChargeQueue = evsChargeQueue;
    }

    void ConstructMappings(const grid_info& grid_i, const order_info& order_i) {
        Orders = GetOrdersMapping(order_i);
        Grids = GetGridsMapping(grid_i);
    }

private:
    void UpdateOrdersState(const EV_info& ev_i, const order_info& order_i) {
        for (i64 evId = 0; evId < ev_i.N_EV; evId++) {
            {
                auto it = EvTargetOrder[evId].begin();
                while (it != EvTargetOrder[evId].end()) {
                    if (!order_i.state[Orders[*it]]) {
                        it++;
                        continue;
                    }
                    EvCarringOrder[evId].insert(*it);
                    EvTargetOrder[evId].erase(it++);
                }
            }

            {
                auto it = EvCarringOrder[evId].begin();
                while (it != EvCarringOrder[evId].end()) {
                    if (Orders.count(*it)) {
                        it++;
                        continue;
                    }
                    EvCarringOrder[evId].erase(it++);
                }
            }
        }

        for (i64 evId = 0; evId < ev_i.N_EV; evId++) {
            for (auto orderId : EvTargetOrder[evId]) {
                AssignedOrders.erase(orderId);
            }
            EvTargetOrder[evId].clear();
        }
    }

    void AssignOrders(const EV_info& ev_i, const order_info& order_i) {
        set<pair<i64, pair<i64, i64>>> ordersQueue;
        for (i64 i = 0; i < order_i.N_order; i++) {
            if (order_i.state[i]) {
                continue;
            }

            i64 orderId = order_i.id[i];
            if (AssignedOrders.count(orderId)) {
                continue;
            }

            for (auto [_, evId] : EvsQueue) {
                if ((i64)EvTargetOrder[evId].size() + (i64)EvCarringOrder[evId].size() >= EV.N_trans_max) {
                    continue;
                }

                i64 cost = EsimateDeliveryCost(ev_i, order_i, Orders, evId, orderId);
                if (cost > T_max) {
                    continue;
                }
                ordersQueue.insert({ cost, {evId, orderId} });
            }
        }

        while (!ordersQueue.empty()) {
            auto it = ordersQueue.begin();
            auto diff = it->first;
            auto [evId, orderId] = it->second;

            if (AssignedOrders.count(orderId)) {
                ordersQueue.erase(it++);
                continue;
            }

            if ((i64)EvTargetOrder[evId].size() + (i64)EvCarringOrder[evId].size() >= EV.N_trans_max) {
                ordersQueue.erase(it++);
                continue;
            }

            i64 cost = EsimateDeliveryCost(ev_i, order_i, Orders, evId, orderId);
            if (diff != cost) {
                ordersQueue.erase(it++);
                ordersQueue.insert({ cost, {evId, orderId} });
                continue;
            }

            EvTargetOrder[evId].insert(orderId);
            AssignedOrders.insert(orderId);
            ordersQueue.erase(it++);
        }
    }

public:
    void ClearChargingPlan(i64 evId) {
        for (auto& [gridId, charges] : ConsumptionByEv[evId]) {
            for (auto [ftm, charge] : charges) {
                ConsumptionByGrid[gridId][ftm] -= charge;
            }
        }
        ConsumptionByEv[evId].clear();
    }


    i64 CalculateSafeCharge(const grid_info& grid_i, i64 gridId, i64 curTm, i64 startTm) {
        // save some time for easy case
        if (true) {
            i64 c = grid_i.y[gridId];

            for (i64 qtm = curTm; qtm < startTm; qtm++) {
                i64 delta = GetPw(grid_i, gridId, qtm, curTm) - ConsumptionByGrid[gridId][qtm];

                i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - c);
                if (excess < 0) {
                    excess = 0;
                }

                i64 buy = -delta - min(grid.V_grid_max, c);
                if (buy < 0) {
                    buy = 0;
                }

                c += delta;
                if (c < 0) {
                    c = 0;
                }
                if (c > grid.C_grid_max) {
                    c = grid.C_grid_max;
                }
            }

            i64 fullCharge = EV.C_EV_max * 4;
            i64 safeCharge = 0;

            for (i64 qtm = startTm; qtm < T_max && fullCharge > 0; qtm++) {
                i64 curCharge = min(fullCharge, EV.V_EV_max);

                i64 delta = GetPw(grid_i, gridId, qtm, curTm) - ConsumptionByGrid[gridId][qtm];

                i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - c);
                if (excess < 0) {
                    excess = 0;
                }

                i64 buy = -delta - min(grid.V_grid_max, c);

                if (buy + curCharge > 0) {
                    curCharge = min(-buy, curCharge);
                    safeCharge += curCharge;
                    fullCharge -= curCharge;
                    break;
                }

                if (buy < 0) {
                    buy = 0;
                }

                safeCharge += curCharge;
                fullCharge -= curCharge;

                c += delta - curCharge;
                if (c < 0) {
                    c = 0;
                }
                if (c > grid.C_grid_max) {
                    c = grid.C_grid_max;
                }
            }

            i64 minCharge = grid.C_grid_max;
            i64 maxCharge = 0;

            auto [excess, buy] = CalculateLosses(grid_i, gridId, curTm, startTm + (safeCharge / EV.V_EV_max + (safeCharge % EV.V_EV_max ? 1 : 0)) - curTm, nullptr, &minCharge, &maxCharge);
            if (!excess && !buy) {
                return min(safeCharge, minCharge);
            }
        }

        if (true) {
            i64 l = 0;
            i64 r = EV.C_EV_max * 4 + 1;

            auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, curTm, startTm - curTm);

            while (r - l > 1) {
                i64 m = (l + r) / 2;

                SimulateCharging(gridId, tm, m, false);
                auto [excess, buy] = CalculateLosses(grid_i, gridId, curTm, startTm - curTm);
                SimulateCharging(gridId, tm, m, true);

                if (buy > etalonBuy) {
                    r = m;
                }
                else {
                    l = m;
                }
            }

            return l;
        }
    }

    i64 CalculateSafeDrop(const grid_info& grid_i, i64 gridId, i64 curTm, i64 startTm) {
        // save some time for easy case
        if (true) {
            i64 c = grid_i.y[gridId];

            for (i64 qtm = curTm; qtm < startTm; qtm++) {
                i64 delta = GetPw(grid_i, gridId, qtm, curTm) - ConsumptionByGrid[gridId][qtm];

                i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - c);
                if (excess < 0) {
                    excess = 0;
                }

                i64 buy = -delta - min(grid.V_grid_max, c);
                if (buy < 0) {
                    buy = 0;
                }

                c += delta;
                if (c < 0) {
                    c = 0;
                }
                if (c > grid.C_grid_max) {
                    c = grid.C_grid_max;
                }
            }

            i64 fullDrop = EV.C_EV_max * 4;
            i64 safeDrop = 0;

            for (i64 qtm = startTm; qtm < T_max && fullDrop > 0; qtm++) {
                i64 curDrop = min(fullDrop, EV.V_EV_max);

                i64 delta = GetPw(grid_i, gridId, qtm, curTm) - ConsumptionByGrid[gridId][qtm];

                i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - c);

                if (excess + curDrop > 0) {
                    curDrop = min(-excess, curDrop);
                    safeDrop += curDrop;
                    fullDrop -= curDrop;
                    break;
                }

                if (excess < 0) {
                    excess = 0;
                }

                i64 buy = -delta - min(grid.V_grid_max, c);

                if (buy < 0) {
                    buy = 0;
                }

                safeDrop += curDrop;
                fullDrop -= curDrop;

                c += delta + curDrop;

                if (c < 0) {
                    c = 0;
                }
                if (c > grid.C_grid_max) {
                    c = grid.C_grid_max;
                }
            }

            i64 minCharge = grid.C_grid_max;
            i64 maxCharge = 0;

            auto [excess, buy] = CalculateLosses(grid_i, gridId, curTm, startTm + (safeDrop / EV.V_EV_max + (safeDrop % EV.V_EV_max ? 1 : 0)) - curTm, nullptr, &minCharge, &maxCharge, true);
            if (!excess && !buy) {
                return min(safeDrop, grid.C_grid_max - maxCharge);
            }
        }

        if (true) {
            i64 l = 0;
            i64 r = EV.C_EV_max * 4 + 1;

            while (r - l > 1) {
                i64 m = (l + r) / 2;

                SimulateCharging(gridId, startTm, m, true);
                auto [excess, buy] = CalculateLosses(grid_i, gridId, curTm, startTm - curTm);
                SimulateCharging(gridId, startTm, m, false);

                auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, curTm, startTm - curTm);
                if (excess > etalonExcess) {
                    r = m;
                }
                else {
                    l = m;
                }
            }

            return l;
        }
    }

public:
    void MatchEvsToGrid(const grid_info& grid_i, const EV_info& ev_i) {
        set<pair<double, pair<i64, i64>>> q;

        for (i64 evId = 0; evId < EV.N_EV; evId++) {
            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                q.insert({ -1e15, {evId, gridId} });
            }
        }

        i64 cnt = 0;

        while (!q.empty()) {
            auto [etVal, t] = *q.begin();
            auto [evId, gridId] = t;
            q.erase(q.begin());

            auto& ev = ev_i.c[evId];
            if (EvChargeGrid[evId] >= 0) {
                continue;
            }

            i64 len = gs.len[ev.u][grid_i.x[gridId]];
            i64 safeCharge = CalculateSafeCharge(grid_i, gridId, tm, tm + len);
            i64 chargeNeeded = EV.C_EV_max - (ev.charge - len * EV.Delta_EV_move);

            if (safeCharge < chargeNeeded && len >= 10) {
                continue;
            }

            chargeNeeded = min(safeCharge, chargeNeeded);

            double val = safeCharge / (len + 1e-8);

            if (-etVal != val) {
                q.insert({ -val, {evId, gridId} });
                continue;
            }

            EvChargeGrid[evId] = grid_i.x[gridId];
            SimulateCharging(gridId, tm + len, chargeNeeded, false, evId);

            cnt += 1;
            if (cnt >= GetEvLimit()) {
                break;
            }
        }
    }

    void HandleEvsQueue(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        for (auto te : EvsQueue) {
            auto evId = te.second;
            auto& ev = ev_i.c[evId];

            // keep moving
            if (ev.u != ev.v) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            if (EvChargeGrid[evId] >= 0) {
                if (ev.u != EvChargeGrid[evId]) {
                    EvTargetGrid[evId] = gs.next[ev.u][EvChargeGrid[evId]];
                    enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                    continue;
                }

                i64 gridId = Grids[EvChargeGrid[evId]];
                if (ConsumptionByEv[evId][gridId].count(tm) && ConsumptionByEv[evId][gridId][tm]) {
                    enqueue(evId, "charge_from_grid " + to_string(ConsumptionByEv[evId][gridId][tm]));
                    continue;
                }

                EvChargeGrid[evId] = -1;
            }

            ClearChargingPlan(evId);

            if (EvTargetGrid[evId] == ev.u) {
                EvTargetGrid[evId] = -1;
            }

            // handle pickup
            if (true) {
                bool pickedUp = false;
                for (i64 orderId : EvTargetOrder[evId]) {
                    if (order_i.w[Orders[orderId]] != ev.u) {
                        continue;
                    }
                    enqueue(evId, "pickup " + to_string(orderId));
                    pickedUp = true;
                    break;
                }

                if (pickedUp) {
                    continue;
                }
            }

            set<i64> dest;
            for (i64 orderId : EvTargetOrder[evId]) {
                dest.insert(order_i.w[Orders[orderId]]);
            }
            for (i64 orderId : EvCarringOrder[evId]) {
                dest.insert(order_i.z[Orders[orderId]]);
            }

            {
                i64 closestDist = (i64)1e18;
                i64 closestGrid = -1;

                for (auto grid : dest) {
                    i64 dist = gs.len[ev.u][grid];
                    if (closestDist > dist) {
                        closestDist = dist;
                        closestGrid = grid;
                    }
                }

                if (closestGrid >= 0) {
                    EvTargetGrid[evId] = gs.next[ev.u][closestGrid];
                }
            }

            if (EvTargetGrid[evId] < 0) {
                EvTargetGrid[evId] = ev.u;
            }

            // true random charge
            if (Grids.count(ev.u) && ((tm < T_last && (grid.DayType == 0 || grid.DayType == 3)) || ev.charge < (T_max - tm) * EV.Delta_EV_move) && ev.charge + EV.V_EV_max <= EV.C_EV_max) {
                i64 gridId = Grids[ev.u];

                i64 safeCharge = min(EV.C_EV_max - ev.charge, CalculateSafeCharge(grid_i, gridId, tm, tm) - EV.V_EV_max);

                i64 minChargeLimit = min(115ll, (T_max - tm)) * EV.Delta_EV_move;
                if (ev.charge < minChargeLimit) {
                    safeCharge = max(safeCharge, minChargeLimit - ev.charge);
                }

                if (safeCharge > EV.V_EV_max) {
                    SimulateCharging(gridId, tm, safeCharge, false, evId);
                    enqueue(evId, "charge_from_grid " + to_string(ConsumptionByEv[evId][gridId][tm]));
                    continue;
                }
            }

            if (EvTargetGrid[evId] == -1 || EvTargetGrid[evId] == ev.u) {
                enqueue(evId, "stay");
                continue;
            }

            enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
        }
    }

public:
    void CalculateIncidents(const grid_info& grid_i, i64 gridId) {
        Incidents[gridId].clear();

        vector<i64> losses(T_max);

        i64 c = grid_i.y[gridId];
        for (i64 qtm = tm; qtm < T_max; qtm++) {
            i64 delta = GetPw(grid_i, gridId, qtm, tm) - ConsumptionByGrid[gridId][qtm];

            i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - c);
            if (excess < 0) {
                excess = 0;
            }

            i64 buy = -delta - min(grid.V_grid_max, c);
            if (buy < 0) {
                buy = 0;
            }

            if (buy > 0 && excess > 0) {
                cerr << "buy > 0 && excess > 0" << endl;
                throw 1;
            }

            losses[qtm] = -buy + excess;

            c += delta;
            if (c < 0) {
                c = 0;
            }
            if (c > grid.C_grid_max) {
                c = grid.C_grid_max;
            }
        }

        i64 sum = 0;
        i64 start = 0;
        for (i64 qtm = tm; qtm < T_max; qtm++) {
            if (!losses[qtm]) {
                if (sum) {
                    Incidents[gridId].push_back({ {start, qtm - 1}, sum });
                }
                sum = 0;
                start = 0;
                continue;
            }

            if (!sum) {
                start = qtm;
            }
            sum += losses[qtm];
        }

        if (sum) {
            Incidents[gridId].push_back({ { start, T_max - 1 }, sum });
        }

        if (true) {
            for (i64 i = 0; i < (i64)Incidents[gridId].size(); i++) {
                auto [t, charge] = Incidents[gridId][i];
                auto [start, end] = t;

                if (start - tm > 200) {
                    Incidents[gridId].resize(i);
                    break;
                }
            }
        }
    }

    void ClearConsumptionByEv(i64 evId) {
        auto it = ConsumptionByEv[evId].begin();
        while (it != ConsumptionByEv[evId].end()) {
            auto tit = it->second.begin();
            while (tit != it->second.end()) {
                if (!tit->second || tit->first < tm) {
                    it->second.erase(tit++);
                }
                else {
                    tit++;
                }
            }
            if (it->second.empty()) {
                ConsumptionByEv[evId].erase(it++);
            }
            else {
                it++;
            }
        }
    }

    void UpdateGridChargePrediction(const grid_info& grid_i, i64 gridId) {
        CalculateIncidents(grid_i, gridId);

        GridLimits[gridId] = { min(-CalculateSafeCharge(grid_i, gridId, tm, tm) + EV.V_EV_max, 0ll), CalculateSafeDrop(grid_i, gridId, tm, tm) };

        GridNextLoss[gridId] = 0;
        for (auto [_, charge] : Incidents[gridId]) {
            if (GridNextLoss[gridId] * charge < 0) {
                break;
            }
            GridNextLoss[gridId] += charge;
        }

        GridChargeNeeded[gridId] = -GridNextLoss[gridId];
        GridChargeNeeded[gridId] = max(GridChargeNeeded[gridId], GridLimits[gridId].first);
        GridChargeNeeded[gridId] = min(GridChargeNeeded[gridId], GridLimits[gridId].second);
    }

    void CalculateLimits(const grid_info& grid_i, const EV_info& ev_i) {
        for (auto [_, evId] : EvsChargeQueue) {
            if (ev_i.c[evId].u != ev_i.c[evId].v) {
                continue;
            }
            ClearChargingPlan(evId);
        }

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            UpdateGridChargePrediction(grid_i, gridId);
        }

        set<pair<pair<double, i64>, pair<i64, i64>>> q;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            if (!GridNextLoss[gridId]) {
                continue;
            }
            for (auto [_, evId] : EvsChargeQueue) {
                auto& ev = ev_i.c[evId];

                if (ev.u != ev.v) {
                    continue;
                }

                if (!ConsumptionByEv[evId].empty()) {
                    continue;
                }

                i64 len = gs.len[ev.u][grid_i.x[gridId]];
                if (ev.charge < len * EV.Delta_EV_move) {
                    continue;
                }

                q.insert({ {-1e9, 0}, {evId, gridId} });
            }
        }

        while (!q.empty()) {
            auto [f, s] = *q.begin();
            q.erase(q.begin());

            auto [etalonVal, _] = f;
            etalonVal *= -1;
            auto [evId, gridId] = s;

            if (etalonVal <= 0) {
                continue;
            }

            if (!GridNextLoss[gridId]) {
                continue;
            }

            auto& ev = ev_i.c[evId];

            if (ev.u != ev.v) {
                continue;
            }

            if (!ConsumptionByEv[evId].empty()) {
                continue;
            }

            double bestVal = 0;
            i64 bestChargeNeeded = 0;
            i64 bestAdditionalGridId = -1;
            i64 bestAdditionalCharge = 0;
            i64 bestSumLen = 0;
            i64 bestDeltaTm = 0;

            for (i64 aGridId = -1; aGridId < grid.N_grid; aGridId++) {
                if (gridId == aGridId) {
                    continue;
                }

                i64 chargeNeeded = GridNextLoss[gridId];

                i64 additionalCharge = 0;
                i64 evCharge = ev.charge;
                i64 deltaTm = 0;

                if (aGridId < 0) {
                    i64 len = gs.len[ev.u][grid_i.x[gridId]];
                    evCharge -= len * EV.Delta_EV_move;
                    if (evCharge < 0) {
                        continue;
                    }

                    deltaTm = len;
                }
                else {
                    i64 aTm = tm;

                    i64 aLen = gs.len[ev.u][grid_i.x[aGridId]];
                    evCharge -= aLen * EV.Delta_EV_move;
                    if (evCharge < 0) {
                        continue;
                    }

                    aTm += aLen;

                    i64 bLen = gs.len[grid_i.x[aGridId]][grid_i.x[gridId]];
                    if (chargeNeeded < 0) {
                        additionalCharge = -GridLimits[aGridId].first;
                        additionalCharge = min(additionalCharge, EV.C_EV_max - evCharge);
                    }
                    else {
                        additionalCharge = GridLimits[aGridId].second;
                        additionalCharge = min(additionalCharge, evCharge - bLen * EV.Delta_EV_move);
                    }

                    if (additionalCharge <= 0) {
                        continue;
                    }

                    aTm += (additionalCharge / EV.V_EV_max) + (additionalCharge % EV.V_EV_max ? 1 : 0);
                    if (chargeNeeded > 0) {
                        additionalCharge *= -1;
                    }
                    evCharge += additionalCharge;

                    evCharge -= bLen * EV.Delta_EV_move;
                    if (evCharge < 0) {
                        continue;
                    }
                    aTm += bLen;

                    deltaTm = aTm - tm;
                }

                i64 chargeAvailable = chargeNeeded < 0 ? evCharge : EV.C_EV_max - evCharge;

                if (chargeNeeded < 0) {
                    chargeNeeded = -min(-chargeNeeded, chargeAvailable);
                }
                else {
                    chargeNeeded = min(chargeNeeded, chargeAvailable);
                }

                auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);
                SimulateCharging(gridId, tm + deltaTm, chargeNeeded, false);
                auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                SimulateCharging(gridId, tm + deltaTm, chargeNeeded, true);

                i64 sumLen = aGridId < 0 ? gs.len[ev.u][grid_i.x[gridId]] : gs.len[ev.u][grid_i.x[aGridId]] + gs.len[grid_i.x[aGridId]][grid_i.x[gridId]];

                double val = -sumLen * EV.Delta_EV_move + gamma * (etalonBuy - buy) + (etalonExcess - excess);
                val *= 10;

                if (grid.DayType == 3) {
                    val += ev.charge - evCharge;
                }

                if (bestVal < val) {
                    bestVal = val;
                    bestChargeNeeded = chargeNeeded;
                    bestAdditionalGridId = aGridId;
                    bestAdditionalCharge = additionalCharge;
                    bestSumLen = sumLen;
                    bestDeltaTm = deltaTm;
                }
            }

            if (bestVal <= 0) {
                continue;
            }

            if (bestVal != etalonVal) {
                q.insert({ {-bestVal, bestSumLen}, {evId, gridId} });
                continue;
            }

            if (bestAdditionalGridId >= 0) {
                SimulateCharging(bestAdditionalGridId, tm + gs.len[ev.u][grid_i.x[bestAdditionalGridId]], bestAdditionalCharge, false, evId);
                UpdateGridChargePrediction(grid_i, bestAdditionalGridId);
            }

            SimulateCharging(gridId, tm + bestDeltaTm, bestChargeNeeded, false, evId);
            UpdateGridChargePrediction(grid_i, gridId);
            ClearConsumptionByEv(evId);
        }
    }

    void UpdateEvChargeGrid(i64 evId) {
        EvChargeGrid[evId] = -1;
        i64 bestTm = T_max + 1;
        for (auto [gridId, t] : ConsumptionByEv[evId]) {
            for (auto [qtm, charge] : t) {
                if (!charge) {
                    continue;
                }
                if (qtm < tm) {
                    continue;
                }
                if (bestTm > qtm) {
                    bestTm = qtm;
                    EvChargeGrid[evId] = grid.x[gridId];
                }
            }
        }
    }

    void HandleCustomEvs(const grid_info& grid_i, const EV_info& ev_i) {
        if (!tm) {
            return;
        }

        bool ok = true;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            CalculateIncidents(grid_i, gridId);
            ok &= Incidents[gridId].empty();
        }

        set<i64> z = { 1, 11, 21 };
        set<i64> w = { 1, 51, 101, 151, 201, 401, 601, 801, 901 };

        if (!ok && ((z.count(tm % 50) && grid.DayType % 2 == 1 && tt >= 0 && WeatherChanges.count(tm / 50)) || w.count(tm))) {
            CalculateLimits(grid_i, ev_i);
        }

        for (auto [_, evId] : EvsChargeQueue) {
            auto& ev = ev_i.c[evId];

            if (ev.u != ev.v) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            UpdateEvChargeGrid(evId);
            if (EvChargeGrid[evId] < 0) {
                continue;
            }

            if (ev.u != EvChargeGrid[evId]) {
                EvTargetGrid[evId] = gs.next[ev.u][EvChargeGrid[evId]];
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            i64 gridId = Grids[EvChargeGrid[evId]];

            if (!ConsumptionByEv[evId][gridId][tm]) {
                continue;
            }

            if (ConsumptionByEv[evId][gridId][tm] > 0) {
                enqueue(evId, "charge_from_grid " + to_string(ConsumptionByEv[evId][gridId][tm]));
                continue;
            }

            if (ConsumptionByEv[evId][gridId][tm] < 0) {
                enqueue(evId, "charge_to_grid " + to_string(-ConsumptionByEv[evId][gridId][tm]));
                continue;
            }
        }
    }

public:
    void DetectWeatherChange(const grid_info& grid_i) {
        i64 div = T_max / grid.N_div;

        i64 cnt = 0;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            i64 patternId = grid.pattern[gridId];
            if (tm % div != 0 && abs(grid.pw_predict[patternId][tm / div] - grid_i.pw_actual[gridId]) > 500) {
                cnt += 1;
            }
        }

        if (cnt == 3) {
            WeatherChanges.insert(tm / div);
        }
    }

public:
    void command(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        if (!tm) {
            cerr << tt << ": " << EvLimit << "/" << EV.N_EV << endl;
        }

        if (tm) {
            DetectWeatherChange(grid_i);
        }

        PreUpdate(grid_i, ev_i, order_i);

        ConstructMappings(grid_i, order_i);

        if (!tm) {
            MatchEvsToGrid(grid_i, ev_i);
        }
        ConstructEvsQueues(ev_i);

        UpdateOrdersState(ev_i, order_i);
        AssignOrders(ev_i, order_i);

        HandleCustomEvs(grid_i, ev_i);
        HandleEvsQueue(grid_i, ev_i, order_i);

        tm += 1;
    }

public:
    pair<double, double> SimulateProccess(map<i64, pair<i64, i64>>& orders, bool& full) {
        grid_info grid_i;
        EV_info ev_i;
        order_info order_i;

        grid_i.N_grid = grid.N_grid;
        grid_i.x.resize(grid_i.N_grid);
        grid_i.y.resize(grid_i.N_grid);
        grid_i.pw_actual.resize(grid_i.N_grid);
        grid_i.pw_excess.resize(grid_i.N_grid);
        grid_i.pw_buy.resize(grid_i.N_grid);

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            grid_i.x[gridId] = grid.x[gridId];
            grid_i.y[gridId] = grid.C_grid_init;
        }

        ev_i.N_EV = EV.N_EV;
        ev_i.c.resize(ev_i.N_EV);

        for (i64 evId = 0; evId < EV.N_EV; evId++) {
            ev_i.c[evId].u = EV.pos[evId];
            ev_i.c[evId].v = EV.pos[evId];
            ev_i.c[evId].dist_from_u = 0;
            ev_i.c[evId].dist_to_v = 0;
            ev_i.c[evId].N_order = 0;
            ev_i.c[evId].charge = EV.C_EV_init;
        }

        order_i.N_order = 0;

        i64 transScore = 0;
        i64 dropped = 0;

        i64 sumPwBuy = 0;
        i64 sumPwExcess = 0;

        i64 maxOrderId = 0;

        for (i64 qtm = 0; qtm <= T_max; qtm++) {
            ConstructMappings(grid_i, order_i);
            UpdateOrdersState(ev_i, order_i);

            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                sumPwBuy += grid_i.pw_buy[gridId];
                sumPwExcess += grid_i.pw_excess[gridId];
            }

            if (qtm == T_max) {
                break;
            }

            if (orders.count(qtm)) {
                order_i.N_order += 1;
                order_i.state.resize(order_i.N_order);
                order_i.time.resize(order_i.N_order);
                order_i.w.resize(order_i.N_order);
                order_i.z.resize(order_i.N_order);
                order_i.id.resize(order_i.N_order);


                order_i.w.back() = orders[qtm].first;
                order_i.z.back() = orders[qtm].second;

                order_i.id.back() = ++maxOrderId;
                order_i.state.back() = 0;
                order_i.time.back() = qtm;
            }

            command(grid_i, ev_i, order_i);
            auto command_list = split_command(dequeue(ev_i));


            vector<i64> consumptionByGrid(grid.N_grid);
            for (i64 evId = 0; evId < EV.N_EV; evId++) {
                auto& ev = ev_i.c[evId];

                stringstream ss(command_list[evId]);

                string cmd;
                ss >> cmd;

                if (cmd == "stay") {
                    continue;
                }

                if (cmd == "pickup") {
                    if (ev.u != ev.v) {
                        cerr << "#" << evId << ": " << cmd << " " << "ev.u != ev.v" << endl;
                        throw 1;
                    }

                    i64 orderId = -1;
                    ss >> orderId;

                    if (!Orders.count(orderId) || order_i.w[Orders[orderId]] != ev.u) {
                        cerr << "#" << evId << ": " << cmd << " " << "there is no such order" << endl;
                        throw 1;
                    }

                    order_i.state[Orders[orderId]] += 1;
                    continue;
                }

                if (cmd == "move") {
                    i64 x = -1;
                    ss >> x;
                    x--;

                    bool canMove = ev.charge >= EV.Delta_EV_move;

                    if (ev.u != ev.v) {
                        if (ev.u != x && ev.v != x) {
                            cerr << "#" << evId << ": " << cmd << " " << "ev.u != x && ev.v != x" << endl;
                            throw 1;
                        }
                        if (canMove) {
                            if (x == ev.u) {
                                ev.dist_from_u -= 1;
                                ev.dist_to_v += 1;
                            }
                            if (x == ev.v) {
                                ev.dist_from_u += 1;
                                ev.dist_to_v -= 1;
                            }
                        }
                    }
                    else {
                        if (!graph.edges[ev.u].count(x)) {
                            cerr << "#" << evId << ": " << cmd << " " << "!graph.edges[ev.u].count(x) " << endl;
                            throw 1;
                        }

                        if (canMove) {
                            ev.v = x;
                            ev.dist_from_u += 1;
                            ev.dist_to_v = graph.edges[ev.u][x] - 1;
                        }
                    }

                    if (!canMove) {
                        continue;
                    }

                    ev.charge -= EV.Delta_EV_move;

                    if (!ev.dist_from_u) {
                        ev.dist_to_v = 0;
                        ev.v = ev.u;
                    }

                    if (!ev.dist_to_v) {
                        ev.dist_from_u = 0;
                        ev.u = ev.v;
                    }

                    if (ev.u == ev.v) {
                        for (auto orderId : EvCarringOrder[evId]) {
                            if (order_i.z[Orders[orderId]] == ev.u) {
                                order_i.state[Orders[orderId]] += 1;
                            }
                        }
                    }
                    continue;
                }

                if (cmd == "charge_from_grid") {
                    i64 v = -1;
                    ss >> v;

                    if (!v) {
                        cerr << "#" << evId << ": " << cmd << " " << "!v" << endl;
                        throw 1;
                    }

                    if (ev.charge + v > EV.C_EV_max) {
                        cerr << "#" << evId << ": " << cmd << " " << "ev.charge + v > EV.C_EV_max" << endl;
                        throw 1;
                    }

                    if (ev.u != ev.v) {
                        cerr << "#" << evId << ": " << cmd << " " << "ev.u != ev.v" << endl;
                        throw 1;
                    }

                    if (!Grids.count(ev.u)) {
                        cerr << "#" << evId << ": " << cmd << " " << "Grids.count(ev.u)" << endl;
                        throw 1;
                    }

                    ev.charge += v;
                    consumptionByGrid[Grids[ev.u]] += v;
                    continue;
                }

                if (cmd == "charge_to_grid") {
                    i64 v = -1;
                    ss >> v;

                    if (!v) {
                        cerr << "#" << evId << ": " << cmd << " " << "!v" << endl;
                        throw 1;
                    }

                    if (ev.charge - v < 0) {
                        cerr << "#" << evId << ": " << cmd << " " << "ev.charge + v > EV.C_EV_max" << endl;
                        throw 1;
                    }

                    if (ev.u != ev.v) {
                        cerr << "#" << evId << ": " << cmd << " " << "ev.u != ev.v" << endl;
                        throw 1;
                    }

                    if (!Grids.count(ev.u)) {
                        cerr << "#" << evId << ": " << cmd << " " << "Grids.count(ev.u)" << endl;
                        throw 1;
                    }

                    ev.charge -= v;
                    consumptionByGrid[Grids[ev.u]] -= v;
                    continue;
                }

                cerr << "#" << evId << ": " << cmd << " " << "invalid cmd" << endl;
                throw 1;
            }

            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                i64& c = grid_i.y[gridId];
                i64 delta = GetPw(grid_i, gridId, qtm, 0) - consumptionByGrid[gridId];
                grid_i.pw_actual[gridId] = GetPw(grid_i, gridId, qtm, 0);

                i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - c);
                if (excess < 0) {
                    excess = 0;
                }
                grid_i.pw_excess[gridId] = excess;

                i64 buy = -delta - min(grid.V_grid_max, c);
                if (buy < 0) {
                    buy = 0;
                }
                grid_i.pw_buy[gridId] = buy;

                c += delta;
                if (c < 0) {
                    c = 0;
                }
                if (c > grid.C_grid_max) {
                    c = grid.C_grid_max;
                }
            }

            if (true) {
                i64 i = 0;
                while (i < order_i.N_order) {
                    if (order_i.state[i] <= 1) {
                        i += 1;
                        continue;
                    }

                    dropped += 1;

                    transScore += T_max - (qtm - order_i.time[i]);

                    swap(order_i.id[i], order_i.id.back());
                    swap(order_i.state[i], order_i.state.back());
                    swap(order_i.time[i], order_i.time.back());
                    swap(order_i.w[i], order_i.w.back());
                    swap(order_i.z[i], order_i.z.back());
                    order_i.N_order -= 1;
                    order_i.id.resize(order_i.N_order);
                    order_i.state.resize(order_i.N_order);
                    order_i.time.resize(order_i.N_order);
                    order_i.w.resize(order_i.N_order);
                    order_i.z.resize(order_i.N_order);
                }
            }
        }

        transScore -= P_trans * order_i.N_order;

        i64 sumEvCharge = 0;
        for (i64 evId = 0; evId < EV.N_EV; evId++) {
            sumEvCharge += ev_i.c[evId].charge;
        }

        i64 sumGridCharge = 0;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            sumGridCharge += grid_i.y[gridId];
        }

        i64 eleScore = (i64)(sumEvCharge + sumGridCharge - gamma * sumPwBuy);

        i64 left = 0;
        i64 picked = 0;
        for (i64 i = 0; i < order_i.N_order; i++) {
            if (!order_i.state[i]) {
                left += 1;
            }
            if (order_i.state[i]) {
                picked += 1;
            }
        }

        cerr << sumPwExcess << "/" << sumPwBuy << "|" << sumEvCharge << "+" << sumGridCharge << "=" << sumEvCharge + sumGridCharge << endl;
        cerr << left << "/" << picked << "/" << dropped << endl;

        cerr << transScore << " " << eleScore << endl;
        cerr << (transScore - S_ref_trans) << " * " << (eleScore - S_ref_ele) << " = " << (transScore - S_ref_trans) * (eleScore - S_ref_ele) << endl;

        full = left + picked == 0;
        if (grid.DayType == 0) {
            full = left + picked < 5;
        }

        return { (double)transScore , (double)eleScore };
    }
};

double CalculateScoreB(vector<pair<double, double>> Ans, const B& prob) {
    Ans.push_back({ (double)prob.S_ref_trans, (double)prob.S_ref_ele });

    for (auto& e : Ans) {
        e.first = max(e.first, (double)prob.S_ref_trans);
        e.second = max(e.second, (double)prob.S_ref_ele);
    }

    std::sort(Ans.begin(), Ans.end());
    std::vector<std::pair<double, double> > a;
    a.push_back(Ans[0]);
    for (size_t i = 1; i < Ans.size(); i++) {
        bool flag = true;
        for (size_t j = i + 1; j < Ans.size(); j++) {
            if (Ans[i].second <= Ans[j].second) {
                flag = false;
                break;
            }
        }
        if (flag) {
            a.push_back(Ans[i]);
        }
    }

    double res = 0;
    for (size_t i = 1; i < a.size(); i++) {
        res += (a[i].first - a[i - 1].first) * (a[i].second - a[0].second);
    }
    return res;
}

void CalculateBestSuite(const vector<pair<pair<double, double>, i64>>& scores, const B& prob, vector<i64>& curSuite, vector<i64>& bestSuite, double& bestScore, const vector<pair<double, double>>& pre) {
    if (pre.size() + curSuite.size() == 5) {
        vector<pair<double, double>> ans = pre;
        for (auto i : curSuite) {
            auto [trans, ele] = scores[i].first;
            i64 delta = 0;
            if (prob.grid.DayType == 3) {
                delta = 300'000;
            }
            ans.push_back({ trans, ele + delta });
        }

        auto curScore = CalculateScoreB(ans, prob);
        if (bestScore < curScore) {
            bestScore = curScore;
            for (i64 j = 0; j < (i64)curSuite.size(); j++) {
                i64 i = curSuite[j];
                bestSuite[pre.size() + j] = scores[i].second;
            }
        }

        return;
    }

    if (prob.grid.DayType != 2) {
        for (i64 c = scores.size() - 1; c >= 0; c--) {
            curSuite.push_back(c);
            CalculateBestSuite(scores, prob, curSuite, bestSuite, bestScore, pre);
            curSuite.pop_back();
        }
    }
    else {
        for (i64 c = 0; c < (i64)scores.size(); c++) {
            curSuite.push_back(c);
            CalculateBestSuite(scores, prob, curSuite, bestSuite, bestScore, pre);
            curSuite.pop_back();
        }
    }
}


int main() {
    setbuf(stderr, nullptr);
    bool isA = false;
    bool dump = false;

    i64 N_solution = 1;
    if (!isA) {
        cin >> N_solution;
    }

    B prob(cin, isA);

    if (dump) {
        stringstream ss;
        prob.dump(ss, isA);

        ofstream ofs("common.dump");
        ofs << ss.str();
        ofs.close();
    }

    TStrategy* str = nullptr;
    graph_summary gs(prob.graph, prob.grid);
    grid_info grid_i(prob.grid.N_grid);
    EV_info ev_i(prob.EV.N_EV);
    order_info order_i;

    string command_per_turn;
    vector<pair<double, double>> scores; scores.reserve(N_solution);
    cerr << "DayType: " << prob.grid.DayType << endl;
    cerr << prob.EV.N_EV << "/" << prob.grid.N_grid << endl;

    double predictedBestScore = 0.0;
    vector<pair<pair<double, double>, i64>> predictedScores;



    if (!isA) {
        BestSuite.resize(5);

        map<i64, pair<i64, i64>> orders;
        if (true) {
            mt19937 rng((unsigned int)(chrono::steady_clock::now().time_since_epoch().count()));
            uniform_int_distribution<i64> dist(0ll, std::numeric_limits<i64>::max());


            const vector<i64> z = { 1, 1, 1, 0, 1, 1, 0, 1, 1, 0 };
            for (i64 qtm = 0; qtm < prob.T_last; qtm++) {
                if (!z[qtm % 10]) {
                    continue;
                }

                i64 w = dist(rng) % prob.graph.V;
                i64 z = w;

                while (true) {
                    z = dist(rng) % prob.graph.V;
                    if (w != z) {
                        break;
                    }
                }

                orders[qtm] = { w, z };
            }
        }

        i64 nmax = prob.EV.N_EV;
        for (i64 i = 1; i < nmax; i++) {
            try {
                str = new TStrategy(prob, gs, -i, isA, scores, { 0, 0 });
                str->initialize();

                bool full = false;
                auto e = str->SimulateProccess(orders, full);
                if (full) {
                    nmax = min(nmax, i + 4);
                }
                predictedScores.push_back({ e, i });
            }
            catch (...) {
                predictedScores.clear();
                break;
            }
        }

        stringstream ss;
        for (auto& e : predictedScores) {
            auto [trans, ele] = e.first;
            ss << (i64)trans << " " << (i64)ele << endl;
        }

        if (dump) {
            ofstream ofs("scores.dump");
            ofs << ss.str();
            ofs.close();
        }
    }

    if (!isA && predictedScores.empty()) {
        BestSuite = { 3, 15, 11, 6, 9 };
    }

    for (i64 n = 0; n < N_solution; ++n) {
        if (!predictedScores.empty() && !isA) {
            vector<i64> curSuite;
            predictedBestScore = 0.0;
            CalculateBestSuite(predictedScores, prob, curSuite, BestSuite, predictedBestScore, scores);

            if (false && prob.grid.DayType == 1) {
                double bestScore = 0;
                for (i64 j = n; j < (i64)BestSuite.size(); j++) {
                    for (auto [t, i] : predictedScores) {
                        auto [trans, ele] = t;
                        if (BestSuite[j] != i) {
                            continue;
                        }
                        scores.push_back(t);
                        auto score = CalculateScoreB(scores, prob);
                        if (bestScore < score) {
                            bestScore = score;
                            BestSuite[n] = i;
                        }
                        scores.pop_back();
                        break;
                    }
                }
            }
        }

        pair<double, double> predictedScore = { 0.0, 0.0 };
        for (auto [t, i] : predictedScores) {
            if ((i64)BestSuite.size() < n || BestSuite[n] != i) {
                continue;
            }
            predictedScore = t;
            break;
        }


        stringstream ss;
        str = new TStrategy(prob, gs, n, isA, scores, predictedScore);
        str->initialize();

        i64 pwExcess = 0;
        i64 pwBuy = 0;

        set<i64> picked;
        set<i64> dropped;
        set<i64> left;

        for (i64 t = 0; t < prob.T_max; ++t) {
            grid_i.load(cin);
            ev_i.load(cin, isA);

            if (!isA) {
                order_i.load(cin);
            }

            if (dump) {
                grid_i.dump(ss);
                ev_i.dump(ss, isA);
                if (!isA) {
                    order_i.dump(ss);
                }
            }

            str->command(grid_i, ev_i, order_i);
            command_per_turn = str->dequeue(ev_i);
            auto command_list = split_command(command_per_turn);
            cout << command_per_turn << flush;

            for (i64 i = 0; i < grid_i.N_grid; i++) {
                pwExcess += grid_i.pw_excess[i];
                pwBuy += grid_i.pw_buy[i];
            }

            {
                set<i64> currentlyPicked;
                for (int i = 0; i < order_i.N_order; i++) {
                    if (!order_i.state[i]) {
                        left.insert(order_i.id[i]);
                        continue;
                    }
                    left.erase(order_i.id[i]);
                    currentlyPicked.insert(order_i.id[i]);
                }

                for (auto id : picked) {
                    if (currentlyPicked.count(id)) {
                        continue;
                    }
                    dropped.insert(id);
                }

                picked = currentlyPicked;
            }
        }


        grid_i.load(cin);
        ev_i.load(cin, isA);
        if (!isA) {
            order_i.load(cin);
        }

        if (dump) {
            grid_i.dump(ss);
            ev_i.dump(ss, isA);
            if (!isA) {
                order_i.dump(ss);
            }
        }

        i64 sumEvCharge = 0;
        for (i64 evId = 0; evId < prob.EV.N_EV; evId++) {
            sumEvCharge += ev_i.c[evId].charge;
        }

        i64 sumGridCharge = 0;
        for (i64 gridId = 0; gridId < prob.grid.N_grid; gridId++) {
            sumGridCharge += grid_i.y[gridId];
        }

        cerr << pwExcess << "/" << pwBuy << "|" << sumEvCharge << "+" << sumGridCharge << "=" << sumEvCharge + sumGridCharge << endl;

        if (!isA) {
            cerr << left.size() << "/" << picked.size() << "/" << dropped.size() << endl;

            double S_trans = 0.0;
            double S_ele = 0.0;
            cin >> S_trans >> S_ele;

            cerr << S_trans << " " << S_ele << endl;
            cerr << (S_trans - prob.S_ref_trans) << " * " << (S_ele - prob.S_ref_ele) << " = " << (S_trans - prob.S_ref_trans) * (S_ele - prob.S_ref_ele) << endl;

            scores.push_back({ S_trans, S_ele });

            if (dump) {
                if (!isA) {
                    ss << S_trans << " " << S_ele << endl;
                }
            }
        }
        else {
            double score = 0.0;
            cin >> score;

            cerr << score << endl;
        }


        if (dump) {
            ofstream ofs;
            ofs.open("sol_" + to_string(n) + ".dump");
            ofs << ss.str();
            ofs.close();
        }

    }

    if (!isA) {
        for (auto [trans, ele] : scores) {
            cerr << (i64)trans << " " << (i64)ele << endl;
        }

        auto score = CalculateScoreB(scores, prob);
        cerr << predictedBestScore << "-" << score << "=" << predictedBestScore - score << endl;
    }

    return 0;
}