#include <iostream>
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

using namespace std;

typedef long long int i64;

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
    i64 DayType;
    i64 N_div, N_pattern, sigma_ele, Delta_event;
    double p_event;
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
    i64 N_EV, C_EV_init, C_EV_max, V_EV_max, N_trans_max, Delta_EV_move;
    std::vector<i64> pos;
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
    i64 N_grid;
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
    i64 N_EV;
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

    //order_info() = default;

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

struct TStrategy : public strategy {
public:
    bool isA = false;

public:
    map<i64, i64> Orders;
    map<i64, i64> Grids;

public:
    set<pair<i64, i64>, greater<pair<i64, i64>>> EvsQueue;
    set<pair<i64, i64>, greater<pair<i64, i64>>> EvsChargeQueue;

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
    vector<i64> PreviousGridCharge;

    i64 SavedExcess = 0;
    i64 WasteExcess = 0;

    vector<i64> PathCnt;

    vector<i64> GridSumGenerated;
    vector<i64> GridSumPredicted;
    vector<i64> GridSumExcess;
    vector<i64> GridSumBuy;

    vector<i64> GridSumDonated;

    vector<i64> EvFullCharge;

    i64 PwActualSum = 0;
    i64 PwActualMax = 0;

    i64 MaxExcess = 0;

public:
    vector<vector<pair<pair<i64, i64>, i64>>> Incidents;

public:
    i64 Succ = 0;
    i64 Miss = 0;

public:
    TStrategy(const B& b, const graph_summary& gs, i64 qt, bool _isA) : strategy(b, gs), tt(qt), isA(_isA) {
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
        EvFullCharge = vector<i64>(EV.N_EV, -1);

        GridSumGenerated = vector<i64>(grid.N_grid);
        GridSumExcess = vector<i64>(grid.N_grid);
        GridSumBuy = vector<i64>(grid.N_grid);
        GridSumDonated = vector<i64>(grid.N_grid);

        Incidents = vector<vector<pair<pair<i64, i64>, i64>>>(grid.N_grid);
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
        if (curTm && tm / div == (curTm - 1) / div && abs(grid.pw_predict[patternId][tm / div] - grid_i.pw_actual[gridId]) > 100) {
            //cerr << tm << " " << grid.pw_predict[patternId][tm / div] << "/" << grid_i.pw_actual[gridId] << endl;
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

    void TestGridsPwPredict(const grid_info& grid_i) {
        if (!tm) {
            return;
        }

        stringstream ss;
        for (i64 i = 0; i < grid_i.N_grid; i++) {
            i64 d = grid_i.pw_actual[i] - GetPw(grid_i, i, tm - 1, tm);
            ss << d << " ";
        }
        cerr << ss.str() << endl;
    }

    void TestGridsPwBalance(const grid_info& grid_i) {
        if (!tm) {
            return;
        }

        for (i64 gridId = 0; gridId < grid_i.N_grid; gridId++) {
            i64 delta = grid_i.pw_actual[gridId] - ConsumptionByGrid[gridId][tm - 1];

            i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - PreviousGridCharge[gridId]);
            if (excess < 0) {
                excess = 0;
            }

            i64 buy = -delta - min(grid.V_grid_max, PreviousGridCharge[gridId]);
            if (buy < 0) {
                buy = 0;
            }

            //if (excess != grid_i.pw_excess[gridId]) {
            //    cerr << tm << " excess != grid_i.pw_excess[gridId]" << endl;
            //    throw 1;
            //}

            //if (buy != grid_i.pw_buy[gridId]) {
            //    cerr << gridId << " buy != grid_i.pw_buy[gridId]" << endl;
            //    throw 1;
            //}
        }
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
    void UpdateChargeHistory(const grid_info& grid_i, const EV_info& ev_i) {
        i64 pwActualCur = 0;
        i64 pwExcessCur = 0;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            pwActualCur += grid_i.pw_actual[gridId];
            pwExcessCur += grid_i.pw_excess[gridId];

            GridSumGenerated[gridId] += grid_i.pw_actual[gridId];
            GridSumExcess[gridId] += grid_i.pw_excess[gridId];
            GridSumBuy[gridId] += grid_i.pw_buy[gridId];
        }
        PwActualSum += pwActualCur;

        if (MaxExcess < pwExcessCur) {
            MaxExcess = pwExcessCur;
        }

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            pwActualCur += grid_i.y[gridId];
        }

        for (i64 evId = 0; evId < EV.N_EV; evId++) {
            pwActualCur += ev_i.c[evId].charge;
        }
        PwActualMax = max(PwActualMax, pwActualCur);
    }

    void PreUpdate(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        PreviousGridCharge = grid_i.y;
        UpdateChargeHistory(grid_i, ev_i);
    }

private:
    void OuputCurrentChargeStat(const grid_info& grid_i, const EV_info& ev_i) {
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            cerr << "grid #" << gridId << ": " << grid_i.y[gridId] << "/" << grid_i.pw_excess[gridId] << endl;
        }
        for (i64 evId = 0; evId < EV.N_EV; evId++) {
            cerr << "ev #" << evId << ": " << ev_i.c[evId].charge << endl;
        }
    }

    void OutputExcessStat(const grid_info& grid_i, const EV_info& ev_i) {
        i64 gridCapacity = 0;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            gridCapacity += grid.C_grid_max - grid_i.y[gridId];
        }

        i64 evCapacity = 0;
        for (auto [_, evId] : EvsChargeQueue) {
            evCapacity += EV.C_EV_max - ev_i.c[evId].charge;
        }

        i64 excess = 0;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            excess += grid_i.pw_excess[gridId];
        }

        cerr << tm << ": " << excess << "<" << gridCapacity << "+" << evCapacity << endl;
    }


    void PreStat(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        OuputCurrentChargeStat(grid_i, ev_i);
        OutputExcessStat(grid_i, ev_i);
    }

private:
    void CalculatePathCnt() {
        PathCnt = vector<i64>(graph.V);

        for (i64 u = 0; u < graph.V; u++) {
            for (i64 v = 0; v < graph.V; v++) {
                if (u == v) {
                    continue;
                }

                i64 c = u;
                while (c != v) {
                    PathCnt[c] += 1;
                    c = gs.next[c][v];
                }
                PathCnt[v] += 1;
            }
        }
    }

    void PredictGeneratedSum(const grid_info& grid_i) {
        GridSumPredicted = vector<i64>(grid.N_grid);
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            for (i64 qtm = 0; qtm < T_max; qtm++) {
                GridSumPredicted[gridId] += GetPw(grid_i, gridId, qtm, tm);
            }
        }
    }

    void FirstUpdate(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        CalculatePathCnt();
        PredictGeneratedSum(grid_i);
    }


private:
    void OutputGridDistances(const grid_info& grid_i) {
        stringstream ss;
        for (i64 i = 0; i < grid_i.N_grid; i++) {
            for (i64 j = 0; j < grid_i.N_grid; j++) {
                ss << setw(4) << gs.len[grid_i.x[i]][grid_i.x[j]] << " ";
            }
            ss << endl;
        }
        cerr << ss.str();
    }

    void OutputGridPatternStat(const grid_info& grid_i) {
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            for (i64 i = 0; i < grid.N_grid; i++) {
                cerr << grid.pw_predict[grid.pattern[gridId]][i] << " ";
            }
            cerr << endl;
        }
    }

    void OutputExcessiveGrid(const grid_info& grid_i) {
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            i64 charge = 0;
            auto [excess, buy] = CalculateLosses(grid_i, gridId, 0, 0, &charge);
            i64 realExcess = excess - (grid.C_grid_max - charge);
            cerr << "gridId #" << gridId << ": " << excess << "/" << charge << "/" << realExcess << endl;
        }
    }

    void FirstStat(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        cerr << tt << ": " << EvLimit << "/" << EV.N_EV << endl;

        //OutputGridDistances(grid_i);
        //OutputGridPatternStat(grid_i);
        //OutputExcessiveGrid(grid_i);
    }

private:
    void OutputChargeSummary(const grid_info& grid_i, const EV_info& ev_i) {
        /*
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            auto gridCharge = grid_i.y[gridId];

            i64 maxChargeSpeed = 0;
            for (i64 qtm = 0; qtm < T_max; qtm++) {
                maxChargeSpeed = max(maxChargeSpeed, GetPw(grid_i, gridId, qtm, 0));
            }

            cerr << "grid #" << gridId << ": " << gridCharge << "/" << GridSumGenerated[gridId] << "/" << GridSumExcess[gridId] << "/" << GridSumDonated[gridId] << "/" << maxChargeSpeed << endl;
        }
        */

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            cerr << "grid #" << gridId << ": " << GridSumExcess[gridId] << "/" << GridSumBuy[gridId] << endl;
        }

        /*
        for (i64 evId = 0; evId < EV.N_EV; evId++) {
            auto evCharge = ev_i.c[evId].charge;
            cerr << "ev #" << evId << ": " << evCharge << endl;
        }
        */
    }

    void LastStat(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        // cerr << "PwActualSum/Max = " << PwActualSum << "/" << PwActualMax << endl;
        OutputChargeSummary(grid_i, ev_i);
    }


private:
    i64 GetEvLimit() {
        if (isA) {
            return 0; // 8; // EV.N_EV / 2;
        }

        vector<i64> evLimits = { 3, 15, 11, 6, 9 }; // todo: replace 12 by 11

        /*
        if (evLimits.front() > EV.N_EV) {
            evLimits.front() = 4; // todo: replace 4 by 3
        }
        */

        return evLimits[tt];
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
        EvsChargeQueue.clear();
        EvsQueue.clear();

        for (i64 evId = 0; evId < ev_i.N_EV; ++evId) {
            if (evId >= EvLimit) {
                EvsChargeQueue.insert({ ev_i.c[evId].charge, evId });
                continue;
            }
            EvsQueue.insert({ ev_i.c[evId].charge, evId });
        }
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
                //cerr << safeCharge << "/" << minCharge << endl;
                Succ++;
                return min(safeCharge, minCharge);
            }
            else {
                Miss++;
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
                //cerr << safeDrop << "/" << grid.C_grid_max - maxCharge << endl;
                Succ++;
                return min(safeDrop, grid.C_grid_max - maxCharge);
            }
            else {
                Miss++;
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
                if (excess > etalonExcess) { //  todo may be comapre with etalon
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
    void HandleEvsChargeQueue(const grid_info& grid_i, const EV_info& ev_i) {
        for (auto [_, evId] : EvsChargeQueue) {
            auto& ev = ev_i.c[evId];

            // keep moving
            if (ev.u != ev.v) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            EvTargetGrid[evId] = -1;

            // clear charging plan
            ClearChargingPlan(evId);

            // charge from excessive grid
            if (EvChargeGrid[evId] == ev.u && Grids.count(ev.u) && ev.charge < EV.C_EV_max) {
                i64 gridId = Grids[ev.u];

                i64 safeCharge = min(EV.C_EV_max - ev.charge, CalculateSafeCharge(grid_i, gridId, tm, tm));

                if (safeCharge >= EV.V_EV_max || (tm > T_max * 5 / 10 && safeCharge > 0)) {
                    SimulateCharging(gridId, tm, safeCharge, false, evId);
                    enqueue(evId, "charge_from_grid " + to_string(ConsumptionByEv[evId][gridId][tm]));
                    continue;
                }
            }

            EvChargeGrid[evId] = -1;

            // find closest dropGrid
            for (i64 qq = 0; qq < 1; qq++) {
                i64 bestDropDist = (i64)1e15;
                i64 bestDropGrid = -1;

                for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                    i64 len = gs.len[ev.u][grid_i.x[gridId]];

                    auto [excess, buy] = CalculateLosses(grid_i, gridId, tm, len);
                    if (excess > 0) {
                        continue;
                    }

                    if (bestDropDist > len) {
                        bestDropDist = len;
                        bestDropGrid = gridId;
                    }
                }

                if (bestDropGrid < 0) {
                    break;
                }

                i64 safeDropCharge = CalculateSafeDrop(grid_i, bestDropGrid, tm, tm + bestDropDist);

                i64 maxDropCharge = min(safeDropCharge, ev.charge - bestDropDist * EV.Delta_EV_move);
                i64 maxDropCnt = maxDropCharge / EV.V_EV_max;
                if (maxDropCnt <= 0) {
                    break;
                }

                i64 bestExcessiveDist = (i64)1e15;
                i64 bestExcessiveGrid = -1;

                for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                    i64 len = gs.len[grid_i.x[bestDropGrid]][grid_i.x[gridId]];
                    auto [excess, buy] = CalculateLosses(grid_i, gridId, tm, bestDropDist + maxDropCnt + len);
                    if (excess <= 0) {
                        continue;
                    }
                    if (bestExcessiveDist > len) {
                        bestExcessiveDist = len;
                        bestExcessiveGrid = gridId;
                    }
                }

                if (bestExcessiveGrid < 0) {
                    break;
                }

                i64 dropCharge = min(safeDropCharge, ev.charge - (bestDropDist + bestExcessiveDist + 5) * EV.Delta_EV_move);
                i64 dropCnt = dropCharge / EV.V_EV_max;
                if (dropCnt <= 0) {
                    break;
                }

                EvTargetGrid[evId] = gs.next[ev.u][grid_i.x[bestDropGrid]];
                SimulateCharging(bestDropGrid, tm + bestDropDist, dropCharge, true, evId);

                if (false) {
                    SimulateCharging(bestExcessiveGrid, tm + bestDropDist + dropCnt + (dropCharge % EV.V_EV_max ? 1 : 0) + bestExcessiveDist, EV.C_EV_max - (ev.charge - (bestDropDist + bestExcessiveDist) * EV.Delta_EV_move - dropCharge), false, evId);
                }
            }

            if (tm > T_max * 9 / 10) {
                EvTargetGrid[evId] = -1;
            }

            if (ev.u == EvTargetGrid[evId]) {
                //cerr << "charge_to_grid " + to_string(-ConsumptionByEv[evId][grids[ev.u]][tm]) << endl;
                enqueue(evId, "charge_to_grid " + to_string(-ConsumptionByEv[evId][Grids[ev.u]][tm]));
                GridSumDonated[Grids[ev.u]] += -ConsumptionByEv[evId][Grids[ev.u]][tm];
                continue;
            }

            if (EvTargetGrid[evId] < 0) {
                if (ev.charge + EV.V_EV_max > EV.C_EV_max) {
                    continue;
                }
                i64 bestExcessiveDist = (i64)1e15;
                i64 bestExcessiveGrid = -1;

                for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                    i64 len = gs.len[ev.u][grid_i.x[gridId]];
                    auto [excess, buy] = CalculateLosses(grid_i, gridId, tm, len);
                    if (excess <= 0) {
                        continue;
                    }
                    if (bestExcessiveDist > len) {
                        bestExcessiveDist = len;
                        bestExcessiveGrid = gridId;
                    }
                }

                if (bestExcessiveGrid < 0) {
                    continue;
                }

                if (ev.charge < bestExcessiveDist * EV.Delta_EV_move) {
                    continue;
                }

                EvChargeGrid[evId] = grid_i.x[bestExcessiveGrid];
                EvTargetGrid[evId] = gs.next[ev.u][grid_i.x[bestExcessiveGrid]];

                if (false) {
                    SimulateCharging(bestExcessiveGrid, tm + bestExcessiveDist, EV.C_EV_max - (ev.charge - bestExcessiveDist * EV.Delta_EV_move), false, evId);
                }
            }

            if (EvTargetGrid[evId] == ev.u) {
                continue;
            }

            enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
        }
    }

    void HandleEvsCharge(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        for (auto [_, evId] : EvsChargeQueue) {
            auto& ev = ev_i.c[evId];

            ClearChargingPlan(evId);
            EvChargeGrid[evId] = -1;

            if (ev.u != ev.v) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            EvTargetGrid[evId] = -1;
        }

        set<pair<i64, pair<pair<i64, i64>, pair<i64, i64>>>> chargingQueue;

        vector<i64> eps(grid.N_grid);
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            for (i64 i = 0; i < grid.N_grid; i++) {
                eps[gridId] = max(eps[gridId], gs.len[grid_i.x[gridId]][grid_i.x[i]]);
            }
        }

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);

            for (auto [_, evId] : EvsChargeQueue) {
                auto& ev = ev_i.c[evId];

                i64 dest = ev.u == ev.v ? ev.u : EvTargetGrid[evId];
                i64 shift = ev.u == dest ? ev.dist_from_u : ev.dist_to_v;

                i64 len = gs.len[dest][grid_i.x[gridId]] + shift;

                // todo: recheck it
                if (ev.charge < len * EV.Delta_EV_move) {
                    continue;
                }

                i64 safeCharge = CalculateSafeCharge(grid_i, gridId, tm, tm + len);
                i64 safeDrop = CalculateSafeDrop(grid_i, gridId, tm, tm + len);

                safeCharge = min(safeCharge, EV.C_EV_max - (ev.charge - len * EV.Delta_EV_move));
                safeDrop = min(safeDrop, (ev.charge - len * EV.Delta_EV_move) - eps[gridId]);

                // todo: recheck it
                if (safeDrop < 0) {
                    continue;
                }

                i64 charge = (safeCharge - safeDrop) / 2;

                if (abs(charge) < EV.V_EV_max) {
                    continue;
                }

                i64 value = -len * EV.Delta_EV_move;

                if (charge > 0) {
                    SimulateCharging(gridId, tm + len, charge, false);
                    auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                    SimulateCharging(gridId, tm + len, charge, true);
                    value += (etalonBuy - buy) * 2 + (etalonExcess - excess);
                }
                else {
                    SimulateCharging(gridId, tm + len, -charge, true);
                    auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                    SimulateCharging(gridId, tm + len, -charge, false);
                    value += (etalonBuy - buy) * 2 + (etalonExcess - excess);
                }

                if (value < 0) {
                    continue;
                }

                chargingQueue.insert({ value, {{len, charge}, {evId, gridId}} });
            }
        }


        while (!chargingQueue.empty()) {
            auto it = chargingQueue.begin();
            i64 queueValue = it->first;
            auto [queueLen, queueCharge] = it->second.first;
            auto [evId, gridId] = it->second.second;

            chargingQueue.erase(it);

            if (EvChargeGrid[evId] >= 0) {
                continue;
            }

            // make this a function
            auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);

            auto& ev = ev_i.c[evId];

            i64 dest = ev.u == ev.v ? ev.u : EvTargetGrid[evId];
            i64 shift = ev.u == dest ? ev.dist_from_u : ev.dist_to_v;

            i64 len = gs.len[dest][grid_i.x[gridId]] + shift;

            // always true
            if (ev.charge < len * EV.Delta_EV_move) {
                continue;
            }

            i64 safeCharge = CalculateSafeCharge(grid_i, gridId, tm, tm + len);
            i64 safeDrop = CalculateSafeDrop(grid_i, gridId, tm, tm + len);

            safeCharge = min(safeCharge, EV.C_EV_max - (ev.charge - len * EV.Delta_EV_move));
            safeDrop = min(safeDrop, (ev.charge - len * EV.Delta_EV_move) - eps[gridId]);

            // always true
            if (safeDrop < 0) {
                continue;
            }

            i64 charge = (safeCharge - safeDrop) / 2;

            if (abs(charge) < EV.V_EV_max) {
                continue;
            }

            i64 value = -len * EV.Delta_EV_move;

            if (charge > 0) {
                SimulateCharging(gridId, tm + len, charge, false);
                auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                SimulateCharging(gridId, tm + len, charge, true);
                value += (etalonBuy - buy) * 2 + (etalonExcess - excess);
            }
            else {
                SimulateCharging(gridId, tm + len, -charge, true);
                auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                SimulateCharging(gridId, tm + len, -charge, false);
                value += (etalonBuy - buy) * 2 + (etalonExcess - excess);
            }

            if (queueValue != value || queueCharge != charge) {
                chargingQueue.insert({ value, {{len, charge}, {evId, gridId}} });
                continue;
            }

            if (charge > 0) {
                SimulateCharging(gridId, tm + len, charge, false, evId);
            }
            else {
                SimulateCharging(gridId, tm + len, -charge, true, evId);
            }

            EvChargeGrid[evId] = grid_i.x[gridId];

            if (ev.u != ev.v) {
                continue;
            }

            if (ev.u != EvChargeGrid[evId]) {
                EvTargetGrid[evId] = gs.next[ev.u][grid_i.x[gridId]];
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            if (ConsumptionByEv[evId][gridId][tm] > 0) {
                enqueue(evId, "charge_from_grid " + to_string(ConsumptionByEv[evId][gridId][tm]));
                continue;
            }
            else if (ConsumptionByEv[evId][gridId][tm] < 0) {
                enqueue(evId, "charge_to_grid " + to_string(-ConsumptionByEv[evId][gridId][tm]));
                continue;
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

            ClearChargingPlan(evId);

            if (EvTargetGrid[evId] == ev.u) {
                EvTargetGrid[evId] = -1;
            }

            // todo: important
            //if (ev.charge < EV.Delta_EV_move && !grids.count(ev.u)) {
            //    cerr << "low charge " << evId << endl;
            //}

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
            if (Grids.count(ev.u) && ev.charge < (T_max - tm) * EV.Delta_EV_move && ev.charge + EV.V_EV_max <= EV.C_EV_max) {
                i64 chargeNeeded = min(EV.C_EV_max, (T_max - tm) * EV.Delta_EV_move);
                i64 gridId = Grids[ev.u];

                /*
                i64 l = 0;
                i64 r = chargeNeeded + 1; // todo : fix this!

                auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);

                while (r - l > 1) {
                    i64 m = (l + r) / 2;

                    SimulateCharging(gridId, tm, m, false, evId);
                    auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                    SimulateCharging(gridId, tm, m, true, evId);

                    //cerr << buy << ">" << etalonBuy << endl;
                    i64 val = (etalonExcess - excess) + (etalonBuy - buy) * 2;
                    if (buy > etalonBuy) {
                        //if (val < 0) {
                        r = m;
                    }
                    else {
                        l = m;
                    }
                }

                i64 safeCharge = l;

                i64 fastSafeCharge = min(chargeNeeded, CalculateSafeCharge(grid_i, gridId, tm, tm));
                if (abs(safeCharge - fastSafeCharge) > 0) {
                    cerr << safeCharge - fastSafeCharge << "=" << safeCharge << "-" << fastSafeCharge << endl;

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

                        c += delta;
                        if (c < 0) {
                            c = 0;
                        }
                        if (c > grid.C_grid_max) {
                            c = grid.C_grid_max;
                        }

                        cerr << qtm << ": " << c << "/" << excess << "/" << buy << "/" << ConsumptionByGrid[gridId][qtm] << "/" << delta << endl;
                    }

                }

                */

                i64 safeCharge = min(EV.C_EV_max - ev.charge, CalculateSafeCharge(grid_i, gridId, tm, tm));
                //i64 safeCharge = min(chargeNeeded, CalculateSafeCharge(grid_i, gridId, tm, tm));
                if (safeCharge > EV.V_EV_max) {
                    SimulateCharging(gridId, tm, safeCharge, false, evId);
                    enqueue(evId, "charge_from_grid " + to_string(ConsumptionByEv[evId][gridId][tm]));
                    continue;
                }
            }

            /*
            // excessive charge
            if (Grids.count(ev.u)) {
                i64 gridId = Grids[ev.u];
                i64 chargeNeeded = EV.C_EV_max - ev.charge;

                i64 l = 0;
                i64 r = chargeNeeded + 1; // todo : fix this!

                auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);

                while (r - l > 1) {
                    i64 m = (l + r) / 2;

                    SimulateCharging(gridId, tm, m, false, evId);
                    auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                    SimulateCharging(gridId, tm, m, true, evId);

                    if (!excess) {
                        r = m;
                    }
                    else {
                        l = m;
                    }
                }

                i64 safeCharge = l;

                if (safeCharge > EV.V_EV_max) {
                    SimulateCharging(gridId, tm, safeCharge, false, evId);
                    enqueue(evId, "charge_from_grid " + to_string(ConsumptionByEv[evId][gridId][tm]));
                    continue;
                }
            }
            */

            if (EvTargetGrid[evId] == -1 || EvTargetGrid[evId] == ev.u) {
                enqueue(evId, "stay");
                continue;
            }

            enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
        }
    }

public:
    void OutputExpectedLosses(const grid_info& grid_i) {
        i64 sumExcess = 0;
        i64 sumBuy = 0;

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
            sumExcess += excess;
            sumBuy += buy;
        }

        cerr << tm << ": " << sumExcess << "/" << sumBuy << endl;
    }

    void OutputPredictedGridState(const grid_info& grid_i) {
        if (tm % 100 != 1) {
            return;
        }
        i64 sumExcess = 0;
        i64 sumBuy = 0;

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {

            map<i64, i64> allExcess;
            map<i64, i64> allBuy;

            vector<i64> losses(T_max);

            i64 c = grid_i.y[gridId];
            for (i64 qtm = tm; qtm < T_max; qtm++) {
                i64 delta = GetPw(grid_i, gridId, qtm, tm) - ConsumptionByGrid[gridId][qtm];

                i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - c);
                if (excess < 0) {
                    excess = 0;
                }

                if (excess > 0) {
                    allExcess[qtm] = excess;
                }

                i64 buy = -delta - min(grid.V_grid_max, c);
                if (buy < 0) {
                    buy = 0;
                }

                if (buy > 0) {
                    allBuy[qtm] = buy;
                }

                if (buy > 0 && excess > 0) {
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

            auto safeCharge = CalculateSafeCharge(grid_i, gridId, tm, tm);
            auto safeDrop = CalculateSafeDrop(grid_i, gridId, tm, tm);

            vector<pair<pair<i64, i64>, i64>> q; // [start; stop] -> charge

            i64 sum = 0;
            i64 start = 0;
            for (i64 qtm = tm; qtm < T_max; qtm++) {
                if (!losses[qtm]) {
                    if (sum) {
                        q.push_back({ {start, qtm - 1}, sum });
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
                q.push_back({ { start, T_max - 1 }, sum });
            }

            if (q.empty()) {
                cerr << "grid #" << gridId << " is fine" << endl;
                continue;
            }

            cerr << "grid #" << gridId << ": ";
            cerr << -safeCharge << "/" << safeDrop << ": ";
            for (auto [t, sum] : q) {
                auto [start, stop] = t;
                if (sum < 0) {
                    sumBuy += -sum;
                }
                else {
                    sumExcess += sum;
                }
                cerr << " [" << start << ", " << stop << "]" << " -> " << sum << ";";
            }
            cerr << endl;

            /*
            if (allExcess.empty() && allBuy.empty()) {
                auto safeCharge = CalculateSafeCharge(grid_i, gridId, tm, tm);
                auto safeDrop = CalculateSafeDrop(grid_i, gridId, tm, tm);

                cerr << "grid #" << gridId << ": " << -safeCharge << "/" << safeDrop << endl;
                continue;
            }

            auto outputLosses = [&gridId](map<i64, i64>& q, i64 m) {
                auto it = q.begin();
                cerr << "grid #" << gridId << "(" << it->first << ")" << ":";
                i64 sum = 0;
                while (it != q.end()) {
                    if (it != q.begin()) {
                        auto pit = it;
                        pit--;
                        if (pit->first + 1 != it->first) {
                            break;
                        }
                    }
                    sum += it->second * m;
                    cerr << " " << it->second * m;
                    it++;
                }
                cerr << " = " << sum;
                cerr << endl;
                return 0;
            };

            if (!allExcess.empty() && !allBuy.empty()) {
                if (allExcess.begin()->first < allBuy.begin()->first) {
                    outputLosses(allExcess, 1);
                }
                else if (allExcess.begin()->first > allBuy.begin()->first) {
                    outputLosses(allBuy, -1);
                }
                else {
                    // excess and buy on the same time is impossible
                    throw 1;
                }
            }
            else if (!allExcess.empty()) {
                outputLosses(allExcess, 1);
            }
            else if (!allBuy.empty()) {
                outputLosses(allBuy, -1);
            }
            */
        }

        cerr << sumExcess << "/" << sumBuy << endl;
    }

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
    }

    bool IsBuying(i64 gridId) {
        for (auto [_, charge] : Incidents[gridId]) {
            if (charge > 0) {
                return true;
            }
        }
        return false;
    }

    bool IsExcessive(i64 gridId) {
        for (auto [_, charge] : Incidents[gridId]) {
            if (charge < 0) {
                return true;
            }
        }
        return false;
    }

    void MatchMixedGrids(const grid_info& grid_i, const EV_info& ev_i) {
        set<pair<i64, pair<i64, i64>>> q;

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            if (!IsBuying(gridId) || !IsExcessive(gridId)) {
                continue;
            }

            for (i64 evId = 0; evId < EV.N_EV; evId++) {
                auto& ev = ev_i.c[evId];
                i64 dest = ev.u == ev.v ? ev.u : EvTargetGrid[evId];
                i64 shift = ev.u == EvTargetGrid[evId] ? ev.dist_from_u : ev.dist_to_v;
                q.insert({ gs.len[dest][grid_i.x[gridId]] + shift, {evId, gridId} });
            }
        }

        for (auto [len, t] : q) {
            auto [evId, gridId] = t;
            if (EvChargeGrid[evId] >= 0) {
                continue;
            }

            if (!IsBuying(gridId) || !IsExcessive(gridId)) {
                continue;
            }

            EvChargeGrid[evId] = grid_i.x[gridId];

            auto& ev = ev_i.c[evId];
            i64 evCharge = ev.charge - len * EV.Delta_EV_move;
            if (evCharge < 0) {
                continue;
            }

            i64 gridCharge = grid_i.y[gridId];
            for (i64 qtm = tm; qtm < T_max; qtm++) {
                i64 charge = 0;
                i64 delta = GetPw(grid_i, gridId, qtm, tm) - ConsumptionByGrid[gridId][qtm];

                i64 excess = delta - min(grid.V_grid_max, grid.C_grid_max - gridCharge);
                if (excess < 0) {
                    excess = 0;
                }
                if (excess > 0 && qtm >= tm + len) {
                    charge = min(excess, min(EV.V_EV_max, EV.C_EV_max - evCharge));
                }

                i64 buy = -delta - min(grid.V_grid_max, gridCharge);
                if (buy < 0) {
                    buy = 0;
                }
                if (buy > 0 && qtm >= tm + len) {
                    charge = max(-buy, -min(EV.V_EV_max, evCharge));
                }

                //if (charge) {
                //    cerr << "#" << gridId << "/" << "#" << evId << ": " << gridCharge << "/" << evCharge << ": " << delta << "/" << charge << endl;
                //}

                evCharge += charge;
                ConsumptionByEv[evId][gridId][qtm] += charge;
                ConsumptionByGrid[gridId][qtm] += charge;
                delta -= charge;

                gridCharge += delta;
                if (gridCharge < 0) {
                    gridCharge = 0;
                }
                if (gridCharge > grid.C_grid_max) {
                    gridCharge = grid.C_grid_max;
                }
            }

            CalculateIncidents(grid_i, gridId);
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

    void CalculateLimits(const grid_info& grid_i, const EV_info& ev_i) {
        vector<pair<i64, i64>> GridLimits(grid.N_grid);
        vector<i64> GridNextLoss(grid.N_grid);
        vector<i64> GridChargeNeeded(grid.N_grid);

        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            GridLimits[gridId] = { -CalculateSafeCharge(grid_i, gridId, tm, tm), CalculateSafeDrop(grid_i, gridId, tm, tm) };

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

        /*
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            cerr << gridId << ": " << "[" << GridLimits[gridId].first << ", " << GridLimits[gridId].second << "] ";
            cerr << GridNextLoss[gridId] << "/" << GridChargeNeeded[gridId] << endl;
        }
        */

        vector<i64> eps(grid.N_grid);
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            for (i64 i = 0; i < grid.N_grid; i++) {
                eps[gridId] = max(eps[gridId], gs.len[grid_i.x[gridId]][grid_i.x[i]]);
            }
        }

        set<pair<pair<i64, i64>, pair<i64, i64>>> q;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            if (!GridNextLoss[gridId]) {
                continue;
            }
            for (auto [_, evId] : EvsChargeQueue) {
                auto& ev = ev_i.c[evId];

                if (ev.u != ev.v) {
                    continue;
                }

                ClearConsumptionByEv(evId);
                if (!ConsumptionByEv[evId].empty()) {
                    continue;
                }

                i64 len = gs.len[ev.u][grid_i.x[gridId]];
                if (ev.charge < len * EV.Delta_EV_move) {
                    continue;
                }

                i64 chargeNeeded = GridNextLoss[gridId];
                if (chargeNeeded < 0) {
                    i64 chargeAvailable = ev.charge - (len + eps[gridId]) * EV.Delta_EV_move;
                    if (chargeAvailable <= 0) {
                        continue;
                    }
                    chargeNeeded = -min(-chargeNeeded, chargeAvailable);
                }
                else {
                    i64 chargeAvailable = EV.C_EV_max - (ev.charge - len * EV.Delta_EV_move);
                    chargeNeeded = min(chargeNeeded, chargeAvailable);
                }

                if (!chargeNeeded) {
                    continue;
                }

                auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);
                SimulateCharging(gridId, tm + len, chargeNeeded, false);
                auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                SimulateCharging(gridId, tm + len, chargeNeeded, true);

                i64 val = -len * EV.Delta_EV_move + 2 * (etalonBuy - buy) + (etalonExcess - excess);
                //if (chargeNeeded < 0) {
                //    cerr << val << endl;
                //}
                q.insert({ {-val, len}, {evId, gridId} });
            }
        }

        while (!q.empty()) {
            auto [f, s] = *q.begin();
            q.erase(q.begin());

            auto [etalonVal, len] = f;
            etalonVal *= -1;
            auto [evId, gridId] = s;
            // cerr << "#" << evId << "->" << "#" << gridId << ": " << etalonVal << "/" << len << endl;

            if (etalonVal <= 0) {
                continue;
            }

            // get rid of copypaste #1
            auto& ev = ev_i.c[evId];

            if (ev.u != ev.v) {
                continue;
            }

            ClearConsumptionByEv(evId);
            if (!ConsumptionByEv[evId].empty()) {
                continue;
            }

            i64 chargeNeeded = GridNextLoss[gridId];
            if (chargeNeeded < 0) {
                i64 chargeAvailable = ev.charge - (len + eps[gridId]) * EV.Delta_EV_move;
                if (chargeAvailable <= 0) {
                    continue;
                }
                chargeNeeded = -min(-chargeNeeded, chargeAvailable);
            }
            else {
                i64 chargeAvailable = EV.C_EV_max - (ev.charge - len * EV.Delta_EV_move);
                chargeNeeded = min(chargeNeeded, chargeAvailable);
            }

            auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);
            SimulateCharging(gridId, tm + len, chargeNeeded, false);
            auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
            SimulateCharging(gridId, tm + len, chargeNeeded, true);

            i64 val = -len * EV.Delta_EV_move + 2 * (etalonBuy - buy) + (etalonExcess - excess);
            if (val != etalonVal) {
                q.insert({ {-val, len}, {evId, gridId} });
                continue;
            }

            //cerr << "#" << evId << "->" << "#" << gridId << ": " << val << "/" << len << "/" << chargeNeeded << endl;

            SimulateCharging(gridId, tm + len, chargeNeeded, false, evId);
            EvChargeGrid[evId] = grid_i.x[gridId];

            // get rid of copypaste #2
            GridLimits[gridId] = { -CalculateSafeCharge(grid_i, gridId, tm, tm), CalculateSafeDrop(grid_i, gridId, tm, tm) };

            CalculateIncidents(grid_i, gridId);
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
    }

    void HandleCustomEvs(const grid_info& grid_i, const EV_info& ev_i) {
        if (!tm) {
            return;
        }

        if (tm % (T_max / grid.N_div) == 1) {
            //cerr << tm << endl;
            //for (i64 evId = 0; evId < EV.N_EV; evId++) {
            //    cerr << "#" << evId << ": " << ev_i.c[evId].charge << endl;
            //}
            //for (auto [_, evId] : EvsChargeQueue) {
            //    ClearChargingPlan(evId);
            //}

            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                CalculateIncidents(grid_i, gridId);

                /*
                cerr << "grid #" << gridId << ":";
                for (auto [t, charge] : Incidents[gridId]) {
                    auto [start, end] = t;
                    cerr << " " << "[" << start << ", " << end << "]" << " -> " << charge << ";";
                }
                cerr << endl;
                */
            }

            CalculateLimits(grid_i, ev_i);

            //MatchMixedGrids(grid_i, ev_i);
        }

        for (auto [_, evId] : EvsChargeQueue) {
            auto& ev = ev_i.c[evId];
            if (EvChargeGrid[evId] < 0) {
                continue;
            }

            if (ev.u != ev.v) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
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

            //cerr << tm << ": " << "#" << gridId << "/" << "#" << evId << ": " << grid_i.y[gridId] << "/" << ev.charge << "/" << ConsumptionByEv[evId][gridId][tm] << endl;

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
    void command(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        if (!tm) {
            FirstStat(grid_i, ev_i, order_i);
            //cerr << gs.diameter << endl;
        }

        TestGridsPwBalance(grid_i);
        PreUpdate(grid_i, ev_i, order_i);

        ConstructMappings(grid_i, order_i);
        ConstructEvsQueues(ev_i);

        UpdateOrdersState(ev_i, order_i);
        AssignOrders(ev_i, order_i);

        //OutputPredictedGridState(grid_i);

        //HandleEvsChargeQueue(grid_i, ev_i);
        //HandleEvsCharge(grid_i, ev_i, order_i);
        HandleCustomEvs(grid_i, ev_i);
        HandleEvsQueue(grid_i, ev_i, order_i);

        tm += 1;

        if (tm == T_max) {
            LastStat(grid_i, ev_i, order_i);
        }
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

int main() {
    setbuf(stderr, nullptr);
    bool isA = false;
    bool dump = true;

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
    for (i64 n = 0; n < N_solution; ++n) {
        stringstream ss;
        str = new TStrategy(prob, gs, n, isA);
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

        if (dump) {
            ofstream ofs;
            ofs.open("sol_" + to_string(n) + ".dump");
            ofs << ss.str();
            ofs.close();
        }

        i64 sumEvCharge = 0;
        for (i64 evId = 0; evId < prob.EV.N_EV; evId++) {
            sumEvCharge += ev_i.c[evId].charge;
        }

        i64 sumGridCharge = 0;
        for (i64 gridId = 0; gridId < prob.grid.N_grid; gridId++) {
            sumGridCharge += grid_i.y[gridId];
        }

        //cerr << prob.grid.N_grid << "/" << prob.EV.N_EV << endl;
        cerr << pwExcess << "/" << pwBuy << "|" << sumEvCharge << "+" << sumGridCharge << "=" << sumEvCharge + sumGridCharge << endl;
        cerr << "Succ/Miss = " << str->Succ << "/" << str->Miss << endl;

        if (!isA) {
            cerr << left.size() << "/" << picked.size() << "/" << dropped.size() << endl;

            double S_trans = 0.0;
            double S_ele = 0.0;
            cin >> S_trans >> S_ele;

            cerr << S_trans << " " << S_ele << endl;
            cerr << (S_trans - prob.S_ref_trans) << " * " << (S_ele - prob.S_ref_ele) << " = " << (S_trans - prob.S_ref_trans) * (S_ele - prob.S_ref_ele) << endl;
        }
        else {
            double score = 0.0;
            cin >> score;

            cerr << score << endl;
        }
    }
    return 0;
}
