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
    EV_data(std::istream& src) {
        src >> N_EV >> C_EV_init >> C_EV_max >> V_EV_max >> N_trans_max >> Delta_EV_move;
        pos.resize(N_EV);
        for (i64 i = 0; i < N_EV; ++i) {
            src >> pos[i];
            --pos[i];
        }
    }

    void dump(ostream& ss) {
        ss << N_EV << " " << C_EV_init << " " << C_EV_max << " " << V_EV_max << " " << N_trans_max << " " << Delta_EV_move << endl;
        for (i64 i = 0; i < N_EV; ++i) {
            ss << pos[i] << endl;
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

    void dump(ostream& ss) {
        graph.dump(ss);
        grid.dump(ss);
        EV.dump(ss);

        ss << p_const_trans << " " << T_last << endl;
        ss << P_trans << " " << gamma << " " << S_ref_ele << " " << S_ref_trans << endl;
        ss << T_max << endl;
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

    void dump(ostream& ss) {
        ss << charge << endl;
        ss << u << " " << v << " " << dist_from_u << " " << dist_to_v << endl;

        ss << N_adj << endl;
        for (i64 i = 0; i < N_adj; ++i) {
            ss << a[i] << " ";
        }
        ss << endl;

        ss << N_order << endl;
        for (i64 i = 0; i < N_order; ++i) {
            ss << o[i] << " ";
        }
        ss << endl;
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
    void dump(ostream& ss) {
        for (i64 i = 0; i < N_EV; ++i) {
            c[i].dump(ss);
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

    void dump(ostream& ss) {
        ss << N_order << endl;
        for (i64 i = 0; i < N_order; ++i) {
            ss << id[i] << " " << w[i] << " " << z[i] << " " << state[i] << " " << time[i] << endl;
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

    i64 SavedExcess = 0;
    i64 WasteExcess = 0;

    vector<i64> PathCnt;

    vector<i64> GridSumGenerated;
    vector<i64> GridSumPredicted;
    vector<i64> GridSumExcess;
    vector<i64> GridSumDonated;

    vector<i64> EvFullCharge;

    i64 PwActualSum = 0;
    i64 PwActualMax = 0;

    i64 MaxExcess = 0;

    map<i64, i64> RejectedChargingEvs;
public:
    GreedyTransportStrategy(const B& b, const graph_summary& gs, i64 qt) : strategy<B>(b, gs), tt(qt) {
        ConsumptionByGrid = vector<vector<i64>>(b.grid.N_grid, vector<i64>(T_max));
        ConsumptionByEv = vector<map<i64, map<i64, i64>>>(EV.N_EV);
    }

    void initialize() {
        strategy::initialize();
        AssignedOrders.clear();
        EvTargetOrder = vector<set<i64>>(EV.N_EV);
        EvCarringOrder = vector<set<i64>>(EV.N_EV);
        EvTargetGrid = vector<i64>(EV.N_EV, -1);
        EvChargeGrid = vector<i64>(EV.N_EV, -1);
        EvFullCharge = vector<i64>(EV.N_EV, -1);

        GridSumGenerated = vector<i64>(grid.N_grid);
        GridSumExcess = vector<i64>(grid.N_grid);
        GridSumDonated = vector<i64>(grid.N_grid);
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

    /*
    i64 GetPw(const grid_info& grid_i, i64 gridId, i64 tm) {
        i64 patternId = grid.pattern[gridId];
        return grid.pw_predict[patternId][tm / grid.N_div];
    }
    */

    i64 GetPw(const grid_info& grid_i, i64 gridId, i64 tm, i64 curTm) {
        i64 patternId = grid.pattern[gridId];
        i64 div = T_max / grid.N_div;
        if (curTm && tm / div == (curTm - 1) / div && abs(grid.pw_predict[patternId][tm / div] - grid_i.pw_actual[gridId]) > 100) {
            //cerr << tm << " " << grid.pw_predict[patternId][tm / div] << "/" << grid_i.pw_actual[gridId] << endl;
            return grid_i.pw_actual[gridId];
        }
        return grid.pw_predict[patternId][tm / div];
    }


    pair<i64, i64> CalculateLosses(const grid_info& grid_i, i64 gridId, i64 curTm, i64 deltaTm = 0, i64* charge = nullptr) {
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

            if (tm - curTm >= deltaTm) {
                sumExcess += excess;
                sumBuy += buy;
            }
        }

        if (charge) {
            *charge = c;
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

    void command(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        //i64 evLimit = EV.N_EV - (EV.N_EV * tt * 3 / 16);
        //vector<i64> evLimits = {15, 12, 9, 6, 3};
        //vector<i64> evLimits = { 18, 15, 12, 9, 4 };
        //vector<i64> evLimits = { 18, 4, 15, 9, 12 };
        //vector<i64> evLimits = { 12, 18, 15, 9, 4 };
        vector<i64> evLimits = { 18, 15, 12, 6, 9 };

        if (evLimits.front() > EV.N_EV) {
            evLimits.front() = 4;
        }
        //reverse(evLimits.begin(), evLimits.end());

        i64 evLimit = evLimits[tt];

        i64 pwActualCur = 0;
        i64 pwExcessCur = 0;
        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            pwActualCur += grid_i.pw_actual[gridId];
            pwExcessCur += grid_i.pw_excess[gridId];
        }
        PwActualSum += pwActualCur;

        if (MaxExcess < pwExcessCur) {
            MaxExcess = pwExcessCur;
            if (false) {
                cerr << tm << ": " << MaxExcess << endl;
                for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                    cerr << "grid #" << gridId << ": " << grid_i.y[gridId] << "/" << grid_i.pw_excess[gridId] << endl;
                }
                for (i64 evId = 0; evId < EV.N_EV; evId++) {
                    cerr << "ev #" << evId << ": " << ev_i.c[evId].charge << endl;
                }
            }
        }


        for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
            pwActualCur += grid_i.y[gridId];
        }

        for (i64 evId = 0; evId < EV.N_EV; evId++) {
            pwActualCur += ev_i.c[evId].charge;
        }
        PwActualMax = max(PwActualMax, pwActualCur);


        if (!tm) {
            cerr << tt << ": " << evLimit << "/" << EV.N_EV << endl;
        }

        if (false && !tm) {
            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                //vector<i64> q;
                //for (i64 i = 0; i < grid.N_div; i++) {
                //    q.push_back(grid.pw_predict[grid.pattern[gridId]][i]);
                //}

                //vector<i64> ps(q.size());
                //ps.back() = q.back();
                //for (i64 i = grid.N_div - 2; i >= 0; i--) {
                //    ps[i] = ps[i + 1] + q[i];
                //}

                //for (i64 i = 0; i < grid.N_div; i++) {
                //    cerr << ps[i] << " ";
                //}
                //cerr << endl;

                for (i64 i = 0; i < grid.N_grid; i++) {
                    cerr << grid.pw_predict[grid.pattern[gridId]][i] << " ";
                }
                cerr << endl;
            }
        }

        if (!tm) {
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

            GridSumPredicted = vector<i64>(grid.N_grid);
            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                for (i64 qtm = 0; qtm < T_max; qtm++) {
                    GridSumPredicted[gridId] += GetPw(grid_i, gridId, qtm, tm);
                }
            }
        }

        if (tm) {
            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                GridSumGenerated[gridId] += grid_i.pw_actual[gridId]; //  grid_i.pw_actual[gridId] > 0 ? grid_i.pw_actual[gridId] : 0; //  GetPw(grid_i, gridId, tm, 0);
                GridSumExcess[gridId] += grid_i.pw_excess[gridId];
                //if (gridId == 15) {
                //    cerr << grid_i.y[gridId] << " + " << grid_i.pw_actual[gridId] << "/" << grid_i.pw_excess[gridId] << endl;
                //}
            }
        }

        map<i64, i64> orders;
        for (i64 i = 0; i < order_i.N_order; i++) {
            orders[order_i.id[i]] = i;
        }

        map<i64, i64> grids;
        for (i64 i = 0; i < grid_i.N_grid; i++) {
            grids[grid_i.x[i]] = i;
        }


        if (false && !tm) {
            stringstream ss;
            for (i64 i = 0; i < grid_i.N_grid; i++) {
                for (i64 j = 0; j < grid_i.N_grid; j++) {
                    ss << setw(4) << gs.len[grid_i.x[i]][grid_i.x[j]] << " ";
                }
                ss << endl;
            }
            cerr << ss.str();
        }

        set<pair<i64, i64>, greater<pair<i64, i64>>> evsQueue;
        set<pair<i64, i64>, greater<pair<i64, i64>>> evsChargeQueue;
        set<pair<i64, i64>> evsStaticQueue;

        {
            auto it = RejectedChargingEvs.begin();
            while (it != RejectedChargingEvs.end()) {
                if (tm - it->second <= 10 || ev_i.c[it->first].u != ev_i.c[it->first].v) {
                    it++;
                    continue;
                }
                EvTargetGrid[it->first] = -1;
                EvChargeGrid[it->first] = -1;
                RejectedChargingEvs.erase(it++);
            }
        }

        for (i64 evId = 0; evId < ev_i.N_EV; ++evId) {
            //if (evId >= 19) {
            //    evsStaticQueue.insert({ ev_i.c[evId].charge, evId });
            //    continue;
            //}
            if (evId >= evLimit && !RejectedChargingEvs.count(evId) && !ev_i.c[evId].N_order) {
                //evsStaticQueue.insert({ ev_i.c[evId].charge, evId });
                evsChargeQueue.insert({ ev_i.c[evId].charge, evId });
                continue;
            }
            evsQueue.insert({ ev_i.c[evId].charge, evId });
        }

        if (false && !tm) {
            map<i64, i64> excessiveGrids;
            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                i64 charge = 0;
                auto [excess, buy] = CalculateLosses(grid_i, gridId, 0, 0, &charge);
                i64 realExcess = excess - (grid.C_grid_max - charge);
                cerr << "gridId #" << gridId << ": " << excess << "/" << charge << "/" << realExcess << endl;
                if (realExcess > 0) {
                    excessiveGrids[gridId] = realExcess;
                }
            }

            set<pair<pair<i64, i64>, pair<i64, i64>>> q;

            for (auto [gridId, excess] : excessiveGrids) {
                for (auto [_, evId] : evsStaticQueue) {
                    i64 len = gs.len[ev_i.c[evId].u][grid_i.x[gridId]];
                    q.insert({ {excess, len}, {evId, gridId} });
                }
            }

            for (auto [key, val] : q) {
                auto [excess, len] = key;
                auto [evId, gridId] = val;

                if (EvChargeGrid[evId] >= 0) {
                    continue;
                }
                if (!excessiveGrids.count(gridId)) {
                    continue;
                }

                EvChargeGrid[evId] = grid_i.x[gridId];
                excessiveGrids.erase(gridId);
            }
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

            //if (evLimit < 10 && gs.len[order_i.w[i]][order_i.z[i]] > gs.diameter / 2) {
            //    continue;
            //}

            //for (i64 evId = 0; evId < ev_i.N_EV; evId++) {
            for (auto [_, evId] : evsQueue) {
                if ((i64)EvTargetOrder[evId].size() + (i64)EvCarringOrder[evId].size() >= EV.N_trans_max) {
                    continue;
                }

                if (RejectedChargingEvs.count(evId) && tm - RejectedChargingEvs[evId] > 10) {
                    continue;
                }

                //if (evLimit < 10 && gs.len[ev_i.c[evId].u][order_i.w[i]] > gs.diameter / 2) {
                //    continue;
                //}

                //if (evId >= evLimit) {
                //    continue;
                //}

                /*
                i64 shift = ev_i.c[evId].u == ev_i.c[evId].v ? 0 : (ev_i.c[evId].u == EvTargetGrid[evId] ? ev_i.c[evId].dist_from_u : ev_i.c[evId].dist_to_v);
                i64 dest = ev_i.c[evId].u == ev_i.c[evId].v ? ev_i.c[evId].u : EvTargetGrid[evId];
                if (evLimit >= 10) {
                    ordersQueue.insert({ shift + (i64)gs.len[dest][order_i.w[i]], {evId, orderId} });
                }
                else {
                    ordersQueue.insert({ shift + (i64)gs.len[dest][order_i.w[i]] + (i64)sqrt(gs.len[order_i.w[i]][order_i.z[i]]), {evId, orderId} });
                }
                */
                i64 cost = EsimateDeliveryCost(ev_i, order_i, orders, evId, orderId);
                if (cost > T_max) {
                    continue;
                }
                ordersQueue.insert({ cost, {evId, orderId} });
            }
        }

        //if (T_max - tm < 100) {
        //    ordersQueue.clear();
        //}

        if (false) {
            for (auto& o : ordersQueue) {
                i64 evId = o.second.first;
                i64 orderId = o.second.second;

                if (AssignedOrders.count(orderId)) {
                    continue;
                }

                if ((i64)EvTargetOrder[evId].size() + (i64)EvCarringOrder[evId].size() >= EV.N_trans_max) {
                    continue;
                }

                if (evId >= evLimit && !RejectedChargingEvs.count(evId)) {
                    continue;
                }

                if (evId >= evLimit && (i64)EvTargetOrder[evId].size() + (i64)EvCarringOrder[evId].size() >= 1) {
                    continue;
                }

                EvTargetOrder[evId].insert(orderId);
                AssignedOrders.insert(orderId);
            }
        }

        if (true) {
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

                i64 cost = EsimateDeliveryCost(ev_i, order_i, orders, evId, orderId);
                if (diff != cost) {
                    ordersQueue.erase(it++);
                    ordersQueue.insert({ cost, {evId, orderId} });
                    continue;
                }

                //cerr << tm << ": " << evId << "/" << orderId << "/" << cost << endl;
                EvTargetOrder[evId].insert(orderId);
                AssignedOrders.insert(orderId);
                ordersQueue.erase(it++);
            }
        }


        if (false) {
            //for (i64 evId = 0; evId < EV.N_EV; evId++) {
            for (auto& [_, evId] : evsQueue) {
                for (auto& [gridId, charges] : ConsumptionByEv[evId]) {
                    for (auto [ftm, charge] : charges) {
                        ConsumptionByGrid[gridId][ftm] -= charge;
                    }
                }
                ConsumptionByEv[evId].clear();
            }
        }


        for (auto [_, evId] : evsChargeQueue) {
            auto& ev = ev_i.c[evId];

            //if (ev.N_order) {
            //    cerr << "ev #" << evId << " is carring an order" << endl;
            //    continue;
            //}

            // keep moving
            if (ev.u != ev.v) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }
            EvTargetGrid[evId] = -1;

            // clear charging plan
            if (true) {
                for (auto& [gridId, charges] : ConsumptionByEv[evId]) {
                    for (auto [ftm, charge] : charges) {
                        ConsumptionByGrid[gridId][ftm] -= charge;
                    }
                }
                ConsumptionByEv[evId].clear();
            }

            // charge from excessive grid
            if (EvChargeGrid[evId] == ev.u && grids.count(ev.u) && ev.charge < EV.C_EV_max) {
                i64 gridId = grids[ev.u];
                i64 chargeNeeded = EV.C_EV_max - ev.charge;

                i64 l = 0;
                i64 r = chargeNeeded + 1;

                auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);

                while (r - l > 1) {
                    i64 m = (l + r) / 2;

                    SimulateCharging(gridId, tm, m, false, evId);
                    auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
                    SimulateCharging(gridId, tm, m, true, evId);

                    if (buy > etalonBuy) {
                        r = m;
                    }
                    else {
                        l = m;
                    }
                }

                i64 safeCharge = l;

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
                    if (excess > 0) { //&& grid_i.y[gridId] < grid.C_grid_max * 3 / 4) {
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

                i64 safeDropCharge = 0;
                {
                    i64 len = gs.len[ev.u][grid_i.x[bestDropGrid]];

                    i64 l = 0;
                    i64 r = grid.C_grid_max + 1;

                    while (r - l > 1) {
                        i64 m = (l + r) / 2;

                        SimulateCharging(bestDropGrid, tm + bestDropDist, m, true, evId);
                        auto [excess, buy] = CalculateLosses(grid_i, bestDropGrid, tm, len);
                        SimulateCharging(bestDropGrid, tm + bestDropDist, m, false, evId);

                        if (excess > 0) {
                            r = m;
                        }
                        else {
                            l = m;
                        }
                    }

                    safeDropCharge = l;
                }

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
                enqueue(evId, "charge_to_grid " + to_string(-ConsumptionByEv[evId][grids[ev.u]][tm]));
                GridSumDonated[grids[ev.u]] += -ConsumptionByEv[evId][grids[ev.u]][tm];
                continue;
            }

            if (EvTargetGrid[evId] < 0) {
                if (ev.charge + EV.V_EV_max > EV.C_EV_max) {
                    RejectedChargingEvs.insert({ evId, tm });
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
                    RejectedChargingEvs.insert({ evId, tm });
                    continue;
                }

                if (ev.charge < bestExcessiveDist * EV.Delta_EV_move) {
                    RejectedChargingEvs.insert({ evId, tm });
                    continue;
                }

                EvChargeGrid[evId] = grid_i.x[bestExcessiveGrid];
                EvTargetGrid[evId] = gs.next[ev.u][grid_i.x[bestExcessiveGrid]];

                if (false) {
                    SimulateCharging(bestExcessiveGrid, tm + bestExcessiveDist, EV.C_EV_max - (ev.charge - bestExcessiveDist * EV.Delta_EV_move), false, evId);
                }
            }

            if (EvTargetGrid[evId] == ev.u) {
                RejectedChargingEvs.insert({ evId, tm });
                continue;
            }

            enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
        }

        if (tm > T_max * 0 / 10) {
            RejectedChargingEvs.clear();
        }
        for (auto [evId, qtm] : RejectedChargingEvs) {
            if (qtm == tm) {
                EvTargetGrid[evId] = -1;
                EvChargeGrid[evId] = -1;
                EvTargetOrder[evId].clear();
            }
        }

        for (auto te : evsQueue) {
            auto evId = te.second;
            auto& ev = ev_i.c[evId];

            //if (evsChargeQueue.count({ ev.charge, evId })) {
            //    continue;
            //}

            //if (evId >= evLimit) {
            //    cerr << "ev #" << evId << " will try to deliver an order" << endl;
            //}

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

            // todo: important
            //if (ev.charge < EV.Delta_EV_move && !grids.count(ev.u)) {
            //    cerr << "low charge " << evId << endl;
            //}

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

        //set<pair<i64, i64>> evsRejectedQueue;



        /*
        for (auto [_, evId] : evsStaticQueue) {
            auto& ev = ev_i.c[evId];
            if (EvChargeGrid[evId] < 0) {
                continue;
            }
            // keep moving
            if (ev.u != ev.v) {
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            EvTargetGrid[evId] = -1;
            if (ev.u != EvChargeGrid[evId]) {
                EvTargetGrid[evId] = gs.next[ev.u][EvChargeGrid[evId]];
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            i64 gridId = grids[EvChargeGrid[evId]];

            i64 charge = min(EV.C_EV_max - ev.charge, EV.V_EV_max);

            auto [etalonExcess, etalonBuy] = CalculateLosses(grid_i, gridId, tm);
            SimulateCharging(gridId, tm, charge, false, evId);
            auto [excess, buy] = CalculateLosses(grid_i, gridId, tm);
            SimulateCharging(gridId, tm, charge, true, evId);

            if (excess < etalonExcess) {
                SimulateCharging(gridId, tm, charge, false, evId);
                enqueue(evId, "charge_from_grid " + to_string(ConsumptionByEv[evId][gridId][tm]));
                continue;
            }

            i64 pw = GetPw(grid_i, gridId, tm, tm);
            if (pw < 0) {
                i64 drop = min(-pw, min(EV.V_EV_max, ev.charge));
                if (drop > 0) {
                    SimulateCharging(gridId, tm, drop, true, evId);
                    enqueue(evId, "charge_to_grid " + to_string(-ConsumptionByEv[evId][gridId][tm]));
                    continue;
                }
            }
        }
        */

        tm += 1;

        if (false) {
            i64 gridCapacity = 0;
            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                gridCapacity += grid.C_grid_max - grid_i.y[gridId];
            }
            //for (i64 evId = 0; evId < EV.N_EV; evId++) {
            i64 evCapacity = 0;
            for (auto [_, evId] : evsChargeQueue) {
                evCapacity += EV.C_EV_max - ev_i.c[evId].charge;
            }

            i64 excess = 0;
            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                excess += grid_i.pw_excess[gridId];
            }

            cerr << tm << ": " << excess << "<" << gridCapacity << "+" << evCapacity << endl;

        }

        if (false && tm == T_max) {
            cerr << "PwActualSum/Max = " << PwActualSum << "/" << PwActualMax << endl;
        }

        if (false && tm == T_max) {
            /*
            for (auto [_, evId] : evsChargeQueue) {
                if (EvChargeGrid[evId] < 0) {
                    continue;
                }
                auto evCharge = ev_i.c[evId].charge;
                auto gridCharge = grid_i.y[grids[EvChargeGrid[evId]]];
                cerr << evId << " on #" << grids[EvChargeGrid[evId]] << ": " << evCharge << "/" << gridCharge << " full at " << EvFullCharge[evId] << endl;
            }

            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                auto gridCharge = grid_i.y[gridId];
                double value = (double)(grid.C_grid_init - grid.C_grid_max + GridSumGenerated[gridId]) / PathCnt[grid_i.x[gridId]];
                double predictedValue = (double)(grid.C_grid_init - grid.C_grid_max + GridSumPredicted[gridId]) / PathCnt[grid_i.x[gridId]];

                cerr << "grid #" << gridId << ": " << gridCharge << "/" << PathCnt[grid_i.x[gridId]] << "/" << grid.C_grid_init + GridSumGenerated[gridId] << "/" << GridSumExcess[gridId] << " " << value << "/" << predictedValue << endl;
            }
            */

            for (i64 gridId = 0; gridId < grid.N_grid; gridId++) {
                auto gridCharge = grid_i.y[gridId];

                i64 maxChargeSpeed = 0;
                for (i64 qtm = 0; qtm < T_max; qtm++) {
                    maxChargeSpeed = max(maxChargeSpeed, GetPw(grid_i, gridId, qtm, 0));
                }

                cerr << "grid #" << gridId << ": " << gridCharge << "/" << GridSumGenerated[gridId] << "/" << GridSumExcess[gridId] << "/" << GridSumDonated[gridId] << "/" << maxChargeSpeed << endl;
            }

            for (i64 evId = 0; evId < EV.N_EV; evId++) {
                auto evCharge = ev_i.c[evId].charge;
                cerr << "ev #" << evId << ": " << evCharge << endl;
            }

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
    bool dump = false;
    setbuf(log_dest, nullptr);
    i64 N_solution = 1;
    cin >> N_solution;
    B prob(cin);
    if (dump) {
        stringstream ss;
        prob.dump(ss);

        ofstream ofs("common.dump");
        ofs << ss.str();
        ofs.close();
    }
    std::shared_ptr<strategy<B>> str = nullptr;
    graph_summary gs(prob.graph, prob.grid);
    grid_info grid_i(prob.grid.N_grid);
    EV_info ev_i(prob.EV.N_EV);
    order_info order_i;
    string command_per_turn;
    vector<pair<double, double>> scores; scores.reserve(N_solution);
    cerr << "DayType: " << prob.grid.DayType << endl;
    for (i64 n = 0; n < N_solution; ++n) {
        stringstream ss;
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
            if (dump) {
                grid_i.dump(ss);
                ev_i.dump(ss);
                order_i.dump(ss);
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
        ev_i.load(cin);
        order_i.load(cin);

        if (dump) {
            grid_i.dump(ss);
            ev_i.dump(ss);
            order_i.dump(ss);
        }

        if (dump) {
            ofstream ofs;
            ofs.open("sol_" + to_string(n) + ".dump");
            ofs << ss.str();
            ofs.close();
        }

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
