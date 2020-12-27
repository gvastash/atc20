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

using namespace std;

typedef long long int i64;

FILE* log_dest =
stderr;
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
    grid_data& operator=(const grid_data&) = default;
};
struct EV_data {
    i64 N_EV, C_EV_init, C_EV_max, V_EV_max, N_trans_max, Delta_EV_move;
    std::vector<i64> pos;
    EV_data(std::istream& src) {
        src >> N_EV >> C_EV_init >> C_EV_max >> V_EV_max >> N_trans_max >> Delta_EV_move;
        pos.resize(N_EV);
        for (i64 i = 0; i < N_EV; ++i) {
            src >> pos[i];
            --pos[i];
        }
    }
};
struct A {
    graph_data graph;
    grid_data grid;
    EV_data EV;
    double gamma;
    i64 T_max;
    A(std::istream& src) :graph(src), grid(src), EV(src) {
        for (i64 i = 0; i < grid.N_grid; ++i) {
        }
        src >> gamma;
        src >> T_max;
    }
};
struct B {
    graph_data graph;
    grid_data grid;
    EV_data EV;
    double p_const_trans;
    i64 T_last;
    i64 P_trans;
    double gamma;
    int S_ref_ele, S_ref_trans;
    i64 T_max;
    B(std::istream& src) :graph(src), grid(src), EV(src) {
        for (i64 i = 0; i < grid.N_grid; ++i) {
        }
        src >> p_const_trans >> T_last;
        src >> P_trans >> gamma >> S_ref_ele >> S_ref_trans;
        src >> T_max;
    }
};
struct carinfo {
    i64 charge;
    i64 u, v, dist_from_u, dist_to_v;
    i64 N_adj;
    std::vector<i64> a;
    i64 N_order;
    std::vector<i64> o;
    void load(std::istream& src, [[maybe_unused]] i64 C_EV_max = 25000, [[maybe_unused]] i64 V = 225, [[maybe_unused]] i64 MaxDegree = 5, [[maybe_unused]] i64 N_trans_max = 4, [[maybe_unused]] i64 T_last = 900) {
        src >> charge;
        src >> u >> v >> dist_from_u >> dist_to_v;
        --u, --v;
        src >> N_adj; a.resize(N_adj);
        for (i64 i = 0; i < N_adj; ++i) {
            src >> a[i];
            --a[i];
        }
        src >> N_order; o.resize(N_order);
        for (i64 i = 0; i < N_order; ++i) {
            src >> o[i];
            --o[i];
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
};
std::ostream& operator<<(std::ostream& dest, const grid_info& i) {
    dest << "\tGrid info:\n";
    for (i64 j = 0; j < i.N_grid; ++j)
        dest << "\t\tx: " << i.x[j] << ", y: " << i.y[j] << ", actual: " << i.pw_actual[j] << ", excess: " << i.pw_excess[j] << ", buy: " << i.pw_buy[j] << "\n";
    return dest;
}
struct EV_info {
    i64 N_EV;
    std::vector<carinfo> c;
    EV_info() = default;
    EV_info(i64 N_EV)
        :N_EV(N_EV), c(N_EV) {}
    void load(std::istream& src) {
        for (i64 i = 0; i < N_EV; ++i) {
            c[i].load(src);
        }
    }
};
std::ostream& operator<<(std::ostream& dest, const EV_info& i) {
    dest << "\tEV info:\n";
    for (i64 j = 0; j < i.N_EV; ++j)
        dest << "\t\tcar " << j << "\n";
    return dest;
}
struct order_info {
    i64 N_order;
    std::vector<i64> id, w, z, state, time;
    order_info() = default;
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
};
std::ostream& operator<<(std::ostream& dest, const order_info& i) {
    dest << "\tOrder info: " << i.N_order << " orders left\n";
    for (i64 j = 0; j < i.N_order; ++j)
        dest << "\t\tid: " << i.id[j] << ", departure: " << i.w[j] << ", arrival: " << i.z[j] << ", state: " << i.state[j] << ", ordered at: " << i.time[j] << "\n";
    return dest;
}
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
i64 transit_length(const std::vector<i64>& path, const std::vector<std::vector<i64>>& min_path_len) {
    i64 len = 0;
    for (i64 i = 1; i < (i64)path.size(); ++i)
        len += min_path_len[path[i - 1]][path[i]];
    return len;
}
i64 transit_length(const std::vector<pair<i64, int>>& path, const std::vector<std::vector<i64>>& min_path_len) {
    i64 len = 0;
    for (i64 i = 1; i < (i64)path.size(); ++i)
        len += min_path_len[path[i - 1].first][path[i].first];
    return len;
}
pair<i64, i64> nearest_point(i64 current, const vector<i64>& points, const graph_summary& gs) {
    i64 len = (i64)1e9, nearest_pos = -1, nearest_index = -1;
    for (i64 i = 0; i < (i64)points.size(); ++i)
        if (gs.len[current][points[i]] < len) {
            len = gs.len[current][points[i]];
            nearest_pos = points[i];
            nearest_index = i;
        }
    return { nearest_index, nearest_pos };
}
pair<i64, i64> nearest_nanogrid(i64 current, const graph_summary& gs) {
    return nearest_point(current, gs.nanogrid_pos, gs);
}
string path_string(const vector<pair<i64, int>>& path) {
    string ret;
    for (auto [p, pickup] : path) ret += " -> " + to_string(p + 1) + (pickup != -1 ? "(pickup: " + to_string(pickup) + ")" : "");
    return ret;
}

i64 path_length_test(i64 insert_point, i64 insert_index, const std::vector<i64>& path, const std::vector<std::vector<i64>>& min_path_len) {
    i64 len = insert_index == 0 ? min_path_len[insert_point][path[0]] : 0;
    for (i64 i = 1; i < (i64)path.size(); ++i)
        if (insert_index == i)
            len += min_path_len[path[i - 1]][insert_point] + min_path_len[insert_point][path[i]];
        else
            len += min_path_len[path[i - 1]][path[i]];
    len += insert_index == path.size() ? min_path_len[path.back()][insert_point] : 0;
    return len;
}
struct action : std::list<std::string> {};
struct move_EV : action {
    move_EV(i64 current, i64 goal, const graph_summary& gs) {
        for (i64 cur = current; cur != goal; cur = gs.next[cur][goal]) {
            const i64 next = gs.next[cur][goal];
            for (i64 count = 0; count < gs.len[cur][next]; ++count)
                this->push_back("move " + std::to_string(next + 1));
        }
    }
    move_EV(i64 current, const std::vector<i64>& path, const graph_summary& gs) {
        i64 cur = current;
        for (i64 goal : path)
            for (; cur != goal; cur = gs.next[cur][goal]) {
                const i64 next = gs.next[cur][goal];
                for (i64 count = 0; count < gs.len[cur][next]; ++count)
                    this->push_back("move " + std::to_string(next + 1));
            }
    }
};
auto minimal_matching(const vector<i64>& start, const vector<i64>& goal, const graph_summary& gs) {
    auto minimal_s = start.begin(), minimal_g = goal.begin();
    i64 minimal_len = (i64)1e9;
    for (auto s = start.begin(); s != start.end(); ++s)
        for (auto g = goal.begin(); g != goal.end(); ++g)
            if (gs.len[*s][*g] < minimal_len) {
                minimal_s = s;
                minimal_g = g;
                minimal_len = gs.len[*s][*g];
            }
    return make_pair(minimal_s, minimal_g);
}
template<class... Args>
string strprintf(const char* fmt, const Args & ...args) {
    char buf[65536];
    sprintf(buf, fmt, args...);
    return buf;
}
template<class P>
struct strategy :public P {
    const graph_summary& gs;
    vector<list<string>> command_queue;
    strategy(const P& p, const graph_summary& gs) : P(p), gs(gs),
        command_queue(P::EV.N_EV) {}
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
        //cerr << EV_index << ": " << cmd << endl;
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

struct GreedyTransportStrategy : strategy<B> {
public:
    set<i64> AssignedOrders;

    vector<set<i64>> EvTargetOrder;
    vector<set<i64>> EvCarringOrder;

    vector<i64> EvTargetGrid;
    vector<i64> EvChargeGrid;

    i64 tt = 0;
    i64 tm = 0;

    vector<i64> PreviousGridCharge;

    vector<vector<i64>> ConsumptionByGrid;
    vector<map<i64, map<i64, i64>>> ConsumptionByEv;

public:
    GreedyTransportStrategy(const B& b, const graph_summary& gs, i64 qt) : strategy<B>(b, gs), tt(qt) {
        ConsumptionByGrid = vector<vector<i64>>(b.grid.N_grid, vector<i64>(T_max));
        ConsumptionByEv = vector<map<i64, map<i64, i64>>>(EV.N_EV);
    }

    void initialize() {
        strategy::initialize();
        AssignedOrders.clear();
        EvTargetOrder.resize(EV.N_EV);
        EvCarringOrder.resize(EV.N_EV);
        EvTargetGrid.resize(EV.N_EV, -1);
        EvChargeGrid.resize(EV.N_EV, -1);
    }

    void SimulateCharging(i64 gridId, i64 startTm, i64 chargeNeeded, bool rollback, i64 evId = -1) {
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

    i64 GetPw(const grid_info& grid_i, i64 gridId, i64 tm) {
        i64 patternId = grid.pattern[gridId];
        return grid.pw_predict[patternId][tm / grid.N_div];
    }

    i64 GetPw(const grid_info& grid_i, i64 gridId, i64 tm, i64 curTm) {
        i64 patternId = grid.pattern[gridId];
        i64 div = T_max / grid.N_div;
        if (tm / div == (curTm - 1) / div && abs(grid.pw_predict[patternId][tm / div] - grid_i.pw_actual[gridId]) > 100) {
            //cerr << tm << " " << grid.pw_predict[patternId][tm / div] << "/" << grid_i.pw_actual[gridId] << endl;
            return grid_i.pw_actual[gridId];
        }
        return grid.pw_predict[patternId][tm / div];
    }


    pair<i64, i64> CalculateLosses(const grid_info& grid_i, i64 gridId, i64 curTm) {
        i64 sumExcess = 0;
        i64 sumBuy = 0;

        i64 c = grid_i.y[gridId];

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

            sumExcess += excess;
            sumBuy += buy;
        }

        return { sumExcess, sumBuy };
    }

    i64 CalculateBestGrid(const grid_info& grid_i, const EV_info& ev_i, i64 evId, i64 from, i64 curTm, i64 chargeLeft, i64 buyLimit = 0) {
        auto& ev = ev_i.c[evId];

        double bestValue = -1e18;
        i64 bestGridId = -1;

        for (i64 gridId = 0; gridId < grid_i.N_grid; gridId++) {
            i64 len = gs.len[from][grid_i.x[gridId]];

            if (chargeLeft < len * EV.Delta_EV_move) {
                continue;
            }

            i64 chargeNeeded = (EV.C_EV_max - (chargeLeft - len * EV.Delta_EV_move));
            i64 startTm = curTm + len;

            auto [curExcess, curBuy] = CalculateLosses(grid_i, gridId, curTm);
            SimulateCharging(gridId, startTm, chargeNeeded, /* rollback */ false);

            auto [nxtExcess, nxtBuy] = CalculateLosses(grid_i, gridId, curTm);
            SimulateCharging(gridId, startTm, chargeNeeded, /* rollback */ true);

            i64 excess = nxtExcess - curExcess;
            i64 buy = nxtBuy - curBuy;

            //cerr << gridId << " " << excess << " " << buy << endl;

            //if (buy > buyLimit) {
            //    continue;
            //}

            double value = -len * (i64)EV.Delta_EV_move - buy * gamma - excess;

            if (bestValue < value) {
                bestValue = value;
                bestGridId = gridId;
            }
        }

        //cerr << bestValue << endl;

        return bestGridId >= 0 ? grid_i.x[bestGridId] : -1;
    }

    void TestGridsPwPredict(const grid_info& grid_i) {
        if (!tm) {
            return;
        }

        stringstream ss;
        for (i64 i = 0; i < grid_i.N_grid; i++) {
            i64 d = grid_i.pw_actual[i] - GetPw(grid_i, i, tm - 1);
            ss << d << " ";
        }
        cerr << ss.str() << endl;
    }

    void TestGridsPwBalance(const grid_info & grid_i) {
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

            /*
            if (excess != grid_i.pw_excess[gridId]) {
                cerr << tm << " excess != grid_i.pw_excess[gridId]" << endl;
                throw 1;
            }

            if (buy != grid_i.pw_buy[gridId]) {
                //for (i64 evId = 0; evId < EV.N_EV; evId++) {
                //    if (!ConsumptionByEv[evId].count(gridId)) {
                //        continue;
                //    }
                //    if (!ConsumptionByEv[evId][gridId].count(tm - 1)) {
                //        continue;
                //    }
                //    cerr << evId << " " << ConsumptionByEv[evId][gridId][tm - 1] << endl;
                //}
                //cerr << grid_i.y[gridId] - PreviousGridCharge[gridId] << "/" << grid_i.pw_actual[gridId] << endl;
                //cerr << buy << " - " << grid_i.pw_buy[gridId] << endl;
                //cerr << grid.V_grid_max << "/" << ConsumptionByGrid[gridId][tm - 1] << "/" << grid_i.pw_actual[gridId] << endl;
                //cerr << gridId << " buy != grid_i.pw_buy[gridId]" << endl;
                throw 1;
            }
            */
        }
    }

    void command(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        //i64 evLimit = EV.N_EV - (EV.N_EV * tt * 3 / 16);
        //vector<i64> evLimits = {15, 12, 9, 6, 3};
        //vector<i64> evLimits = { 18, 15, 12, 9, 4 };
        vector<i64> evLimits = { 18, 15, 12, 9, 4 };
        if (evLimits.front() > EV.N_EV) {
            evLimits.front() = 7;
        }
        i64 evLimit = evLimits[tt];


        if (!tm) {
            cerr << tt << ": " << evLimit << "/" << EV.N_EV << endl;
        }

        map<i64, i64> orders;
        for (i64 i = 0; i < order_i.N_order; i++) {
            orders[order_i.id[i]] = i;
        }

        map<i64, i64> grids;
        for (i64 i = 0; i < grid_i.N_grid; i++) {
            grids[grid_i.x[i]] = i;
        }

        set<pair<i64, i64>, greater<pair<i64, i64>>> evsQueue;
        set<pair<i64, i64>, greater<pair<i64, i64>>> evsChargeQueue;

        for (i64 evId = 0; evId < ev_i.N_EV; ++evId) {
            if (evId >= evLimit) {
                evsChargeQueue.insert({ ev_i.c[evId].charge, evId });
                continue;
            }
            evsQueue.insert({ ev_i.c[evId].charge, evId });
        }

        TestGridsPwBalance(grid_i);
        PreviousGridCharge = grid_i.y;

        for (i64 evId = 0; evId < ev_i.N_EV; evId++) {
            {
                auto it = EvTargetOrder[evId].begin();
                while (it != EvTargetOrder[evId].end()) {
                    if (!order_i.state[orders[*it]]) {
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
                    if (orders.count(*it)) {
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

        set<pair<i64, pair<i64, i64>>> ordersQueue;
        for (i64 i = 0; i < order_i.N_order; i++) {
            if (order_i.state[i]) {
                continue;
            }

            i64 orderId = order_i.id[i];
            if (AssignedOrders.count(orderId)) {
                continue;
            }

            //for (i64 evId = 0; evId < ev_i.N_EV; evId++) {
            for (auto [_, evId] : evsQueue) {
                if ((i64)EvTargetOrder[evId].size() + (i64)EvCarringOrder[evId].size() >= EV.N_trans_max) {
                    continue;
                }

                //if (evId >= evLimit) {
                //    continue;
                //}

                i64 shift = ev_i.c[evId].u == ev_i.c[evId].v ? 0 : (ev_i.c[evId].u == EvTargetGrid[evId] ? ev_i.c[evId].dist_from_u : ev_i.c[evId].dist_to_v);
                i64 dest = ev_i.c[evId].u == ev_i.c[evId].v ? ev_i.c[evId].u : EvTargetGrid[evId];
                ordersQueue.insert({ shift + (i64)gs.len[dest][order_i.w[i]], {evId, orderId} });
            }
        }

        //if (T_max - tm < 100) {
        //    ordersQueue.clear();
        //}

        for (auto& o : ordersQueue) {
            i64 evId = o.second.first;
            i64 orderId = o.second.second;

            if (AssignedOrders.count(orderId)) {
                continue;
            }

            if ((i64)EvTargetOrder[evId].size() + (i64)EvCarringOrder[evId].size() >= EV.N_trans_max) {
                continue;
            }

            EvTargetOrder[evId].insert(orderId);
            AssignedOrders.insert(orderId);
        }

        for (auto te : evsQueue) {
            auto evId = te.second;
            auto& ev = ev_i.c[evId];

            if (evId >= evLimit) {
                enqueue(evId, "stay");
                continue;
            }

            if (ev.u != ev.v) {
                /*
                if (EvTargetGrid[evId] != ev.u && EvTargetGrid[evId] != ev.v) {
                    cerr << "Ev " << evId << " on the edge " << ev.u << "-" << ev.v << " cant reach " << EvTargetGrid[evId] << endl;
                    throw 1;
                }

                if (ev.charge < EV.Delta_EV_move) {
                    cerr << "low charge on edge " << evId << endl;
                    throw 1;
                }
                */

                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            // clear charging plan
            if (true) {
                EvChargeGrid[evId] = -1;

                for (auto& [gridId, charges] : ConsumptionByEv[evId]) {
                    for (auto [ftm, charge] : charges) {
                        //cerr << gridId << " at " << ftm << " charge " << charge << endl;
                        ConsumptionByGrid[gridId][ftm] -= charge;
                    }
                }
                ConsumptionByEv[evId].clear();
            }

            if (EvTargetGrid[evId] == ev.u) {
                EvTargetGrid[evId] = -1;
            }

            if (ev.charge < EV.Delta_EV_move && !grids.count(ev.u)) {
                cerr << "low charge " << evId << endl;
            }

            {
                bool pickedUp = false;
                for (i64 orderId : EvTargetOrder[evId]) {
                    if (order_i.w[orders[orderId]] != ev.u) {
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
                dest.insert(order_i.w[orders[orderId]]);
            }
            for (i64 orderId : EvCarringOrder[evId]) {
                dest.insert(order_i.z[orders[orderId]]);
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
            if (grids.count(ev.u) && ev.charge < (T_max - tm) * EV.Delta_EV_move && ev.charge + EV.V_EV_max <= EV.C_EV_max) {
            //if (grids.count(ev.u) && ev.charge + EV.V_EV_max <= EV.C_EV_max) {
                i64 chargeNeeded = min(EV.C_EV_max, (T_max - tm) * EV.Delta_EV_move);
                i64 gridId = grids[ev.u];

                i64 l = 0;
                i64 r = chargeNeeded + 1;

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

                //if (safeCharge > EV.V_EV_max * 10) {
                if (safeCharge > EV.V_EV_max) {
                    SimulateCharging(gridId, tm, safeCharge, false, evId);
                    //cerr << tm << "|" << evId << ": " << safeCharge << "/" << ConsumptionByEv[evId][gridId][tm] << "/"<< ev.charge << endl;
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

        map<i64, i64> excessiveGrids;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
            if (excess > EV.V_EV_max) {
                excessiveGrids[gridId] = excess;
            }
        }

        for (auto [_, evId] : evsChargeQueue) {
            auto& ev = ev_i.c[evId];
            if (ev.u == ev.v && ev.u == EvTargetGrid[evId]) {
                EvTargetGrid[evId] = -1;
                if (EvChargeGrid[evId] >= 0) {
                    EvTargetGrid[evId] = gs.next[ev.u][EvChargeGrid[evId]];
                }
            }
            excessiveGrids.erase(grids[EvChargeGrid[evId]]);
        }

        for (auto [gridId, excess] : excessiveGrids) {
            //cerr << gridId << " " << excess << endl;

            i64 bestDist = (i64)1e16;
            i64 bestEvId = -1;

            for (auto [_, evId] : evsChargeQueue) {
                auto& ev = ev_i.c[evId];

                if (EvTargetGrid[evId] >= 0) {
                    continue;
                }

                if (ev.charge == EV.C_EV_max) {
                    continue;
                }

                i64 dist = min(gs.len[ev.u][grid_i.x[gridId]] + ev.dist_from_u, gs.len[ev.v][grid_i.x[gridId]] + ev.dist_to_v);
                if (bestDist > dist) {
                    bestDist = dist;
                    bestEvId = evId;
                }
            }

            i64 balance = excess - bestDist * EV.Delta_EV_move;
            if (bestEvId >= 0 && balance > 0) {
                //cerr << bestEvId << " will save " << balance << " from grid #" << gridId << endl;
                EvChargeGrid[bestEvId] = grid_i.x[gridId];
                auto& ev = ev_i.c[bestEvId];
                if (ev.u == ev.v) {
                    EvTargetGrid[bestEvId] = gs.next[ev.u][grid_i.x[gridId]];
                }
                else if (gs.len[ev.u][grid_i.x[gridId]] + ev.dist_from_u < gs.len[ev.v][grid_i.x[gridId]] + ev.dist_to_v) {
                    EvTargetGrid[bestEvId] = ev.u;
                }
                else {
                    EvTargetGrid[bestEvId] = ev.v;
                }
            }
        }

        for (auto [_, evId] : evsChargeQueue) {
            auto& ev = ev_i.c[evId];
            if (ev.u != ev.v) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            if (ev.charge == EV.C_EV_max) {
                enqueue(evId, "stay");
                continue;
            }


            if (EvTargetGrid[evId] < 0) {
                enqueue(evId, "stay");
                continue;
            }

            if (ev.u != EvTargetGrid[evId]) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            if (ev.u == EvChargeGrid[evId]) {
                i64 gridId = grids[ev.u];

                i64 c = grid_i.y[gridId];
                i64 delta = GetPw(grid_i, gridId, tm, tm) - ConsumptionByGrid[gridId][tm];

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

                i64 charge = min(excess, min(EV.C_EV_max - ev.charge, EV.V_EV_max));

                if (charge > 0) {
                    //cerr << "save #" << gridId << " by #" << evId << " from excess by charging " << charge << endl;
                    enqueue(evId, "charge_from_grid " + to_string(charge));
                    continue;
                }
            }

            enqueue(evId, "stay");
            continue;
        }

        tm += 1;
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
    setbuf(log_dest, nullptr);
    i64 N_solution = 1;
    cin >> N_solution;
    B prob(cin);
    std::shared_ptr<strategy<B>> str = nullptr;
    graph_summary gs(prob.graph, prob.grid);
    grid_info grid_i(prob.grid.N_grid);
    EV_info ev_i(prob.EV.N_EV);
    order_info order_i;
    string command_per_turn;
    vector<pair<double, double>> scores; scores.reserve(N_solution);
    for (i64 n = 0; n < N_solution; ++n) {
        str.reset(new GreedyTransportStrategy(prob, gs, n));
        str->initialize();

        i64 pwExcess = 0;
        i64 pwBuy = 0;

        set<i64> picked;
        set<i64> dropped;
        set<i64> left;

        for (i64 t = 0; t < prob.T_max; ++t) {
            grid_i.load(cin);
            ev_i.load(cin);
            order_i.load(cin);
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
        ev_i.load(cin);
        order_i.load(cin);
        double S_trans, S_ele;
        cin >> S_trans >> S_ele;


        i64 sumEvCharge = 0;
        for (i64 evId = 0; evId < prob.EV.N_EV; evId++) {
            sumEvCharge += ev_i.c[evId].charge;
        }

        i64 sumGridCharge = 0;
        for (i64 gridId = 0; gridId < prob.grid.N_grid; gridId++) {
            sumGridCharge += grid_i.y[gridId];
        }

        //cerr << prob.grid.N_grid << "/" << prob.EV.N_EV << endl;
        cerr << left.size() << "/" << picked.size() << "/" << dropped.size() << endl;
        cerr << pwExcess << "/" << pwBuy << "|" << sumEvCharge << "+" << sumGridCharge << "=" << sumEvCharge + sumGridCharge << endl;

        cerr << S_trans << " " << S_ele << endl;
        cerr << (S_trans - prob.S_ref_trans) << " * " << (S_ele - prob.S_ref_ele) << " = " << (S_trans - prob.S_ref_trans) * (S_ele - prob.S_ref_ele) << endl;


    }
    return 0;
}
