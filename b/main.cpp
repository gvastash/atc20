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
    constexpr static size_t MaxDegree = 5;
    size_t V, E;
    std::map<size_t, std::map<size_t, size_t>> edges;
    graph_data(std::istream& src) {
        src >> V >> E;
        for (size_t i = 0; i < E; ++i) {
            size_t u, v, d;
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
    size_t DayType;
    size_t N_div, N_pattern, sigma_ele, Delta_event;
    double p_event;
    std::vector<std::vector<int>> pw_predict;
    size_t N_grid, C_grid_init, C_grid_max, V_grid_max;
    std::vector<size_t> x, pattern;
    grid_data() = default;
    grid_data(const grid_data&) = default;
    grid_data(std::istream& src) {
        src >> DayType;
        src >> N_div >> N_pattern >> sigma_ele >> p_event >> Delta_event;
        pw_predict.resize(N_pattern);
        for (size_t i = 0; i < N_pattern; ++i) {
            pw_predict[i].resize(N_div);
            for (size_t j = 0; j < N_div; ++j) {
                src >> pw_predict[i][j];
            }
        }
        src >> N_grid >> C_grid_init >> C_grid_max >> V_grid_max;
        x.resize(N_grid);
        pattern.resize(N_grid);
        for (size_t i = 0; i < N_grid; ++i) {
            src >> x[i] >> pattern[i];
            --x[i];
            --pattern[i];
        }
    }
    grid_data& operator=(const grid_data&) = default;
};
struct EV_data {
    size_t N_EV, C_EV_init, C_EV_max, V_EV_max, N_trans_max, Delta_EV_move;
    std::vector<size_t> pos;
    EV_data(std::istream& src) {
        src >> N_EV >> C_EV_init >> C_EV_max >> V_EV_max >> N_trans_max >> Delta_EV_move;
        pos.resize(N_EV);
        for (size_t i = 0; i < N_EV; ++i) {
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
    size_t T_max;
    A(std::istream& src) :graph(src), grid(src), EV(src) {
        for (size_t i = 0; i < grid.N_grid; ++i) {
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
    size_t T_last;
    size_t P_trans;
    double gamma;
    int S_ref_ele, S_ref_trans;
    size_t T_max;
    B(std::istream& src) :graph(src), grid(src), EV(src) {
        for (size_t i = 0; i < grid.N_grid; ++i) {
        }
        src >> p_const_trans >> T_last;
        src >> P_trans >> gamma >> S_ref_ele >> S_ref_trans;
        src >> T_max;
    }
};
struct carinfo {
    size_t charge;
    size_t u, v, dist_from_u, dist_to_v;
    size_t N_adj;
    std::vector<size_t> a;
    size_t N_order;
    std::vector<size_t> o;
    void load(std::istream& src, [[maybe_unused]] size_t C_EV_max = 25000, [[maybe_unused]] size_t V = 225, [[maybe_unused]] size_t MaxDegree = 5, [[maybe_unused]] size_t N_trans_max = 4, [[maybe_unused]] size_t T_last = 900) {
        src >> charge;
        src >> u >> v >> dist_from_u >> dist_to_v;
        --u, --v;
        src >> N_adj; a.resize(N_adj);
        for (size_t i = 0; i < N_adj; ++i) {
            src >> a[i];
            --a[i];
        }
        src >> N_order; o.resize(N_order);
        for (size_t i = 0; i < N_order; ++i) {
            src >> o[i];
            --o[i];
        }
    }
};
struct grid_info {
    size_t N_grid;
    std::vector<size_t> x, y;
    std::vector<int> pw_actual;
    std::vector<size_t> pw_excess, pw_buy;
    grid_info() = default;
    grid_info(size_t N_grid)
        :N_grid(N_grid), x(N_grid), y(N_grid), pw_actual(N_grid), pw_excess(N_grid), pw_buy(N_grid) {}
    void load(std::istream& src, [[maybe_unused]] size_t V = 225, [[maybe_unused]] size_t C_grid_max = 50000) {
        for (size_t i = 0; i < N_grid; ++i) {
            src >> x[i] >> y[i] >> pw_actual[i] >> pw_excess[i] >> pw_buy[i];
            --x[i];
        }
    }
};
std::ostream& operator<<(std::ostream& dest, const grid_info& i) {
    dest << "\tGrid info:\n";
    for (size_t j = 0; j < i.N_grid; ++j)
        dest << "\t\tx: " << i.x[j] << ", y: " << i.y[j] << ", actual: " << i.pw_actual[j] << ", excess: " << i.pw_excess[j] << ", buy: " << i.pw_buy[j] << "\n";
    return dest;
}
struct EV_info {
    size_t N_EV;
    std::vector<carinfo> c;
    EV_info() = default;
    EV_info(size_t N_EV)
        :N_EV(N_EV), c(N_EV) {}
    void load(std::istream& src) {
        for (size_t i = 0; i < N_EV; ++i) {
            c[i].load(src);
        }
    }
};
std::ostream& operator<<(std::ostream& dest, const EV_info& i) {
    dest << "\tEV info:\n";
    for (size_t j = 0; j < i.N_EV; ++j)
        dest << "\t\tcar " << j << "\n";
    return dest;
}
struct order_info {
    size_t N_order;
    std::vector<size_t> id, w, z, state, time;
    order_info() = default;
    void load(std::istream& src, [[maybe_unused]] size_t V = 225, [[maybe_unused]] size_t T_last = 900) {
        src >> N_order;
        id.resize(N_order);
        w.resize(N_order);
        z.resize(N_order);
        state.resize(N_order);
        time.resize(N_order);
        for (size_t i = 0; i < N_order; ++i) {
            src >> id[i] >> w[i] >> z[i] >> state[i] >> time[i];
            --w[i], --z[i];
        }
    }
};
std::ostream& operator<<(std::ostream& dest, const order_info& i) {
    dest << "\tOrder info: " << i.N_order << " orders left\n";
    for (size_t j = 0; j < i.N_order; ++j)
        dest << "\t\tid: " << i.id[j] << ", departure: " << i.w[j] << ", arrival: " << i.z[j] << ", state: " << i.state[j] << ", ordered at: " << i.time[j] << "\n";
    return dest;
}
struct graph_summary {
    vector<vector<size_t>> len;
    vector<vector<size_t>> next;
    vector<size_t> nanogrid_pos;
    size_t diameter = 0;
    size_t cover_radius = 0;
    graph_summary(const graph_data& graph, const grid_data& grid) :
        len(graph.V, std::vector<size_t>(graph.V, 1e9)),
        next(graph.V, std::vector<size_t>(graph.V)),
        nanogrid_pos(grid.N_grid) {
        const size_t V = graph.V;
        for (size_t i = 0; i < V; ++i)
            len[i][i] = 0;
        for (size_t i = 0; i < V; ++i)
            for (size_t j = 0; j < V; ++j)
                next[i][j] = j;
        for (const auto& [u, u_edges] : graph.edges)
            for (const auto& [v, length] : u_edges) {
                len[u][v] = length;
                len[v][u] = length;
            }
        for (size_t k = 0; k < V; ++k)
            for (size_t i = 0; i < V; ++i)
                for (size_t j = 0; j < V; ++j)
                    if (len[i][j] > len[i][k] + len[k][j]) {
                        len[i][j] = len[i][k] + len[k][j];
                        next[i][j] = next[i][k];
                    }
        nanogrid_pos = grid.x;
        for (size_t i = 0; i < V; ++i)
            for (size_t j = 0; j < V; ++j)
                diameter = max(len[i][j], diameter);
        for (size_t i = 0; i < V; ++i) {
            size_t min_len = 1e9;
            for (size_t j = 0; j < nanogrid_pos.size(); ++j)
                min_len = min(min_len, len[i][j]);
            cover_radius = max(cover_radius, min_len);
        }
    }
};
size_t transit_length(const std::vector<size_t>& path, const std::vector<std::vector<size_t>>& min_path_len) {
    size_t len = 0;
    for (size_t i = 1; i < path.size(); ++i)
        len += min_path_len[path[i - 1]][path[i]];
    return len;
}
size_t transit_length(const std::vector<pair<size_t, int>>& path, const std::vector<std::vector<size_t>>& min_path_len) {
    size_t len = 0;
    for (size_t i = 1; i < path.size(); ++i)
        len += min_path_len[path[i - 1].first][path[i].first];
    return len;
}
pair<size_t, size_t> nearest_point(size_t current, const vector<size_t>& points, const graph_summary& gs) {
    size_t len = 1e9, nearest_pos = -1, nearest_index = -1;
    for (size_t i = 0; i < points.size(); ++i)
        if (gs.len[current][points[i]] < len) {
            len = gs.len[current][points[i]];
            nearest_pos = points[i];
            nearest_index = i;
        }
    return { nearest_index, nearest_pos };
}
pair<size_t, size_t> nearest_nanogrid(size_t current, const graph_summary& gs) {
    return nearest_point(current, gs.nanogrid_pos, gs);
}
string path_string(const vector<pair<size_t, int>>& path) {
    string ret;
    for (auto [p, pickup] : path) ret += " -> " + to_string(p + 1) + (pickup != -1 ? "(pickup: " + to_string(pickup) + ")" : "");
    return ret;
}
vector<pair<size_t, int>> find_transit_path_greedy(size_t current,
    const vector<tuple<size_t, size_t, size_t>>& order,
    const graph_summary& gs) {
    for ([[maybe_unused]] auto [from, to, id] : order) {
    }
    vector<pair<size_t, int>> ret; ret.reserve(2 * order.size());
    vector<size_t> pickup_flag(order.size(), 0);
    vector<size_t> index; index.reserve(order.size());
    vector<size_t> candidate; candidate.reserve(2 * order.size());
    size_t cur = current;
    while (1) {
        for (size_t i = 0; i < order.size(); ++i)
            switch (pickup_flag[i]) {
            case 0:
                candidate.push_back(get<0>(order[i]));
                index.push_back(i);
                break;
            case 1:
                candidate.push_back(get<1>(order[i]));
                index.push_back(i);
                break;
            default:;
            }
        if (candidate.size() == 0) break;
        for ([[maybe_unused]] auto p : candidate);
        auto [i, pos] = nearest_point(cur, candidate, gs);
        ret.emplace_back(pos, pickup_flag[index[i]] == 0 ? get<2>(order[index[i]]) : -1);
        pickup_flag[index[i]] += 1;
        cur = pos;
        candidate.clear();
        index.clear();
    }
    return ret;
}
size_t path_length_test(size_t insert_point, size_t insert_index, const std::vector<size_t>& path, const std::vector<std::vector<size_t>>& min_path_len) {
    size_t len = insert_index == 0 ? min_path_len[insert_point][path[0]] : 0;
    for (size_t i = 1; i < path.size(); ++i)
        if (insert_index == i)
            len += min_path_len[path[i - 1]][insert_point] + min_path_len[insert_point][path[i]];
        else
            len += min_path_len[path[i - 1]][path[i]];
    len += insert_index == path.size() ? min_path_len[path.back()][insert_point] : 0;
    return len;
}
struct action : std::list<std::string> {};
struct move_EV : action {
    move_EV(size_t current, size_t goal, const graph_summary& gs) {
        for (size_t cur = current; cur != goal; cur = gs.next[cur][goal]) {
            const size_t next = gs.next[cur][goal];
            for (size_t count = 0; count < gs.len[cur][next]; ++count)
                this->push_back("move " + std::to_string(next + 1));
        }
    }
    move_EV(size_t current, const std::vector<size_t>& path, const graph_summary& gs) {
        size_t cur = current;
        for (size_t goal : path)
            for (; cur != goal; cur = gs.next[cur][goal]) {
                const size_t next = gs.next[cur][goal];
                for (size_t count = 0; count < gs.len[cur][next]; ++count)
                    this->push_back("move " + std::to_string(next + 1));
            }
    }
};
auto minimal_matching(const vector<size_t>& start, const vector<size_t>& goal, const graph_summary& gs) {
    auto minimal_s = start.begin(), minimal_g = goal.begin();
    size_t minimal_len = 1e9;
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
    bool is_free(size_t EV_index) {
        if (command_queue[EV_index].size() > 0) {
            return false;
        }
        return true;
    }
    string dequeue(const EV_info& ev_i) {
        string ret = "";
        for (size_t i = 0; i < ev_i.N_EV; ++i)
            ret += dequeue(i) + "\n";
        return ret;
    }
    string dequeue(size_t EV_index) {
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
    void enqueue(size_t EV_index, const string& cmd) {
        //cerr << EV_index << ": " << cmd << endl;
        command_queue[EV_index].push_back(cmd);
    }
    void enqueue(size_t EV_index, const string& cmd, size_t repeat) {
        for (size_t i = 0; i < repeat; ++i)
            command_queue[EV_index].push_back(cmd);
    }
    void enqueue(size_t EV_index, list<string>&& cmd_list) {
        command_queue[EV_index].splice(command_queue[EV_index].end(), cmd_list);
    }
};
template<class P>
struct all_stay : strategy<P> {
    all_stay(const P& p, const graph_summary& gs) : strategy<P>(p, gs) {}
    void command(const grid_info&, const EV_info&, const order_info&) {}
};
template<class P>
struct random_walk : strategy<P> {
    using S = strategy<P>;
    std::mt19937_64 engine;
    random_walk(const P& p, const graph_summary& gs) : strategy<P>(p, gs) {}
    void command(const grid_info&, const EV_info& ev_i, const order_info&) {
        for (size_t n = 0; n < ev_i.N_EV; ++n) {
            if (!S::is_free(n)) continue;
            const size_t current = ev_i.c[n].u;
            const size_t safety_energy = S::EV.Delta_EV_move * 50;
            if (auto [_, pos] = nearest_nanogrid(current, S::gs); current != pos) {
                const size_t len_to_charge = S::gs.len[current][pos];
                const int expected_energy = ev_i.c[n].charge - len_to_charge * S::EV.Delta_EV_move;
                if (expected_energy < 0) {
                    S::enqueue(n, "stay", 1000);
                }
                else
                    S::enqueue(n, move_EV(current, pos, S::gs));
                continue;
            }
            else {
                if (ev_i.c[n].charge < safety_energy) {
                    S::enqueue(n, strprintf("charge_from_grid %zu", S::EV.V_EV_max), ceil(1.0 * (safety_energy - ev_i.c[n].charge) / S::EV.V_EV_max));
                    continue;
                }
            }
            uniform_int_distribution<size_t> dice(0, ev_i.c[n].N_adj - 1);
            const size_t goal = dice(engine);
            S::enqueue(n, move_EV(current, ev_i.c[n].a[goal], S::gs));
        }
    }
};

struct less_key {
    const order_info& order_i;
    const graph_summary& gs;
    const int w;

    less_key(const order_info& _order_i, const graph_summary& _gs, int _w) : order_i(_order_i), gs(_gs), w(_w) {

    }

    bool operator() (const int u, const int v)  const {
        long long int from_u = gs.len[w][order_i.w[u]];
        long long int to_u = gs.len[order_i.w[u]][order_i.z[u]];

        long long int from_v = gs.len[w][order_i.w[v]];
        long long int to_v = gs.len[order_i.w[v]][order_i.z[v]];

        double du = (double)from_u + sqrt(to_u);
        double dv = (double)from_v + sqrt(to_v);

        /*
        if (gs.len[w][order_i.w[u]] < gs.len[w][order_i.w[v]]) {
            return true;
        }
        if (gs.len[w][order_i.w[u]] > gs.len[w][order_i.w[v]]) {
            return false;
        }
        */

        if (du < dv) {
            return true;
        }
        if (du > dv) {
            return false;
        }

        return u < v;
    }
};

struct transport_only_0 : strategy<B> {
    std::set<size_t> assigned_order;
    int tt;
    transport_only_0(const B& b, const graph_summary& gs, int qt) :
        strategy<B>(b, gs), tt(qt) {}
    void initialize() {
        strategy::initialize();
        assigned_order.clear();
    }
    void command(const grid_info&, const EV_info& ev_i, const order_info& order_i) {
        for (size_t n = 0; n < ev_i.N_EV; ++n) {
            if (!is_free(n)) continue;
            const size_t current = ev_i.c[n].u;
            const size_t safety_energy = EV.Delta_EV_move * (tt % 2 ? 50 : 100);
            if (auto [_, pos] = nearest_nanogrid(current, gs); current != pos) {
                const size_t len_to_charge = gs.len[current][pos];
                const int expected_energy = ev_i.c[n].charge - len_to_charge * EV.Delta_EV_move;
                if (expected_energy < 0) {
                    enqueue(n, "stay", 1000);
                }
                else
                    enqueue(n, move_EV(current, pos, gs));
                continue;
            }
            else {
                if (ev_i.c[n].charge < safety_energy) {
                    enqueue(n, strprintf("charge_from_grid %zu", EV.V_EV_max), ceil(1.0 * min(EV.C_EV_max, (safety_energy - ev_i.c[n].charge)) / EV.V_EV_max));
                    continue;
                }
            }

            std::vector<size_t> unassigned_order;
            for (size_t i = 0; i < order_i.N_order; ++i)
                if (assigned_order.count(order_i.id[i]) == 0)
                    unassigned_order.push_back(i);

            sort(unassigned_order.begin(), unassigned_order.end(), less_key(order_i, gs, ev_i.c[n].u));

            int i = 0;
            if (i < unassigned_order.size()) {
                size_t count = 0;
                std::vector<tuple<size_t, size_t, size_t>> assign_order;
                while (i < unassigned_order.size() && count++ < EV.N_trans_max) {
                    const size_t order_index = unassigned_order[i];
                    i += 1;

                    const size_t from = order_i.w[order_index];
                    const size_t to = order_i.z[order_index];
                    assign_order.emplace_back(from, to, order_i.id[order_index]);
                    assigned_order.insert(order_i.id[order_index]);
                }
                auto path = find_transit_path_greedy(current, assign_order, gs);
                vector<size_t> transit; transit.reserve(path.size() + 1);
                const size_t expected_transit_length = transit_length(path, gs.len) + gs.len[current][path[0].first];

                size_t radius = tt % 2 ? gs.cover_radius : 0;
                /*
                if (!path.empty()) {
                    for (size_t j = 0; j < gs.nanogrid_pos.size(); ++j) {
                        radius = min(radius, gs.len[0][j]);
                    }
                }
                */

                if (ev_i.c[n].charge < (expected_transit_length + radius) * EV.Delta_EV_move) {
                    enqueue(n, strprintf("charge_from_grid %zu", EV.V_EV_max), (min(EV.C_EV_max, (expected_transit_length + radius) * EV.Delta_EV_move - ev_i.c[n].charge)) / EV.V_EV_max + 1);
                }
                size_t cur = current;
                for (auto [to, pick_up] : path) {
                    enqueue(n, move_EV(cur, to, gs));
                    if (pick_up != -1) enqueue(n, strprintf("pickup %d", pick_up));
                    cur = to;
                }
                continue;
            }
            else {
            }
            continue;
        }
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

public:
    GreedyTransportStrategy(const B& b, const graph_summary& gs, int qt) :
        strategy<B>(b, gs), tt(qt) {}

    void initialize() {
        strategy::initialize();
        AssignedOrders.clear();
        EvTargetOrder.resize(EV.N_EV);
        EvCarringOrder.resize(EV.N_EV);
        EvTargetGrid.resize(EV.N_EV, -1);
        EvChargeGrid.resize(EV.N_EV, -1);
    }

    void command(const grid_info& grid_i, const EV_info& ev_i, const order_info& order_i) {
        tm += 1;

        map<i64, i64> orders;
        for (i64 i = 0; i < order_i.N_order; i++) {
            orders[order_i.id[i]] = i;
        }

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

        set<pair<i64, pair<i64, i64>>> ordersQueue;
        for (i64 i = 0; i < order_i.N_order; i++) {
            if (order_i.state[i]) {
                continue;
            }

            i64 orderId = order_i.id[i];
            if (AssignedOrders.count(orderId)) {
                continue;
            }

            for (i64 evId = 0; evId < ev_i.N_EV; evId++) {
                if (EvTargetOrder[evId].size() + EvCarringOrder[evId].size() >= EV.N_trans_max) {
                    continue;
                }
                //if (EvCarringOrder[evId].size()) {
                //    continue;
                //}

                i64 shift = ev_i.c[evId].u == ev_i.c[evId].v ? 0 : (ev_i.c[evId].u == EvTargetGrid[evId] ? ev_i.c[evId].dist_from_u : ev_i.c[evId].dist_to_v);
                i64 dest = ev_i.c[evId].u == ev_i.c[evId].v ? ev_i.c[evId].u : EvTargetGrid[evId];
                ordersQueue.insert({ shift + (i64)gs.len[dest][order_i.w[i]], {evId, orderId} });
            }
        }

        if (T_max - tm < 100) {
            ordersQueue.clear();
        }

        for (auto& o : ordersQueue) {
            i64 evId = o.second.first;
            i64 orderId = o.second.second;

            if (AssignedOrders.count(orderId)) {
                continue;
            }

            if (EvTargetOrder[evId].size() + EvCarringOrder[evId].size() >= EV.N_trans_max) {
                continue;
            }
            //if (EvCarringOrder[evId].size()) {
            //    continue;
            //}

            //cerr << orderId << " assigned to " << evId << endl;

            EvTargetOrder[evId].insert(orderId);
            AssignedOrders.insert(orderId);
        }


        for (size_t evId = 0; evId < ev_i.N_EV; ++evId) {
            auto& ev = ev_i.c[evId];

            if (ev.u != ev.v) {
                if (EvTargetGrid[evId] != ev.u && EvTargetGrid[evId] != ev.v) {
                    cerr << "Ev " << evId << " on the edge " << ev.u << "-" << ev.v << " cant reach " << EvTargetGrid[evId] << endl;
                    throw 1;
                }
                enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                continue;
            }

            if (EvTargetGrid[evId] == ev.u) {
                EvTargetGrid[evId] = -1;
            }

            if (ev.charge < EV.Delta_EV_move) {
                cerr << "low charge " << evId << endl;
            }

            /*
            if (!EvTargetOrder[evId].empty()) {
                i64 closestDist = 1e18;
                i64 closestOrderId = -1;

                for (auto orderId : EvTargetOrder[evId]) {
                    i64 dist = gs.len[ev.u][order_i.w[orders[orderId]]];

                    if (closestDist > dist) {
                        closestDist = dist;
                        closestOrderId = orderId;
                    }
                }

                EvTargetGrid[evId] = gs.next[ev.u][order_i.w[orders[closestOrderId]]];

                //cerr << ev.u << " -> " << EvTargetGrid[evId] << endl;
                if (ev.u == order_i.w[orders[closestOrderId]]) {
                    enqueue(evId, "pickup " + to_string(closestOrderId));
                    continue;
                }
            }

            if (EvTargetGrid[evId] < 0) {
                if (EvCarringOrder[evId].size() > 0) {
                    EvTargetGrid[evId] = gs.next[ev.u][order_i.z[orders[*EvCarringOrder[evId].begin()]]];
                }
            }
            */

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
                i64 closestDist = 1e18;
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

            bool saveEnergy = true; // tt % 2 == 0;

            // random charge
            {
                i64 charge = 0;
                for (i64 i = 0; i < grid_i.N_grid; i++) {
                    if (grid_i.x[i] != ev.u) {
                        continue;
                    }
                    charge = grid_i.y[i];
                }

                /*
                if (charge < grid.C_grid_max / 4 && ev.charge > EV.C_EV_max / 2) {
                    cerr << "random uncharge" << endl;

                    enqueue(evId, "charge_to_grid " + to_string(EV.V_EV_max));
                    continue;
                }
                */

                if (charge > grid.C_grid_max / 4) {
                    charge = min(charge, (i64)EV.V_EV_max);
                    charge = min(charge, (i64)EV.C_EV_max - (i64)ev.charge);
                    //charge = min(charge, (i64)EV.Delta_EV_move * (i64)(T_max - T_last) - (i64)ev.charge);
                    if (charge > 0) {
                        //cerr << "random charge" << endl;
                        enqueue(evId, "charge_from_grid " + to_string(charge));
                        continue;
                    }
                }
            }

            // check charge
            if (true) {
                i64 len = EvTargetGrid[evId] >= 0 ? gs.len[ev.u][EvTargetGrid[evId]] : -1;

                i64 closestDist = 1e18;
                i64 closestGrid = -1;

                for (auto grid : gs.nanogrid_pos) {

                    bool ok = true;
                    for (i64 i = 0; i < grid_i.N_grid; i++) {
                        if (grid_i.x[i] != grid) {
                            continue;
                        }

                        if (grid_i.y[i] == 0) {
                            ok = false;
                        }
                    }

                    if (!ok && saveEnergy) {
                        continue;
                    }

                    i64 dist = gs.len[EvTargetGrid[evId]][grid];
                    if (closestDist > dist) {
                        closestDist = dist;
                        closestGrid = grid;
                    }
                }

                if ((len + closestDist + 10) * EV.Delta_EV_move >= ev.charge || ev.u == EvTargetGrid[evId] || EvTargetGrid[evId] == -1) {
                    i64 closestDist = 1e18;
                    i64 closestGrid = -1;

                    for (auto grid : gs.nanogrid_pos) {

                        bool ok = true;
                        for (i64 i = 0; i < grid_i.N_grid; i++) {
                            if (grid_i.x[i] != grid) {
                                continue;
                            }

                            if (grid_i.y[i] == 0) {
                                ok = false;
                            }
                        }

                        if (!ok && saveEnergy) {
                            continue;
                        }

                        i64 dist = gs.len[ev.u][grid];
                        if (closestDist > dist) {
                            closestDist = dist;
                            closestGrid = grid;
                        }
                    }

                    EvChargeGrid[evId] = closestGrid;
                }
            }

            // charge
            if (EvChargeGrid[evId] >= 0) {
                if (ev.charge >= EV.C_EV_max) {
                    EvChargeGrid[evId] = -1;
                }
                else {
                    if (ev.u == EvChargeGrid[evId]) {
                        i64 chargeNeeded = min(EV.V_EV_max, EV.C_EV_max - ev.charge);

                        //cerr << evId << " charge_from_grid " + to_string(min(EV.V_EV_max, EV.C_EV_max - ev.charge)) << endl;
                        for (i64 i = 0; i < grid_i.N_grid; i++) {
                            if (grid_i.x[i] != EvChargeGrid[evId]) {
                                continue;
                            }
                            if (saveEnergy) {
                                chargeNeeded = min(chargeNeeded, (i64)grid_i.y[i]);
                            }
                            //cerr << "current grid #" << grid_i.x[i] << " charge = " << grid_i.y[i] << endl;
                        }

                        if (chargeNeeded) {
                            enqueue(evId, "charge_from_grid " + to_string(chargeNeeded));
                        }
                        else {
                            //cerr << "ev #" << evId << " care about charge" << endl;
                            EvChargeGrid[evId] = -1;
                            enqueue(evId, "stay");
                        }

                        continue;
                    }
                    else {
                        EvTargetGrid[evId] = gs.next[ev.u][EvChargeGrid[evId]];
                        enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
                        continue;
                    }
                }
            }

            if (EvTargetGrid[evId] == -1 || EvTargetGrid[evId] == ev.u) {
                enqueue(evId, "stay");
                continue;
            }

            enqueue(evId, "move " + to_string(EvTargetGrid[evId] + 1));
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
enum command_type {
    stay,
    move,
    pickup,
    charge_from_grid,
    charge_to_grid,
    invalid_command
};
struct command {
    command_type type;
    size_t val;
    command(command_type type, size_t val) : type(type), val(val) {}
    string to_str() const {
        switch (type)
        {
        case command_type::stay:
            return strprintf("stay");
        case command_type::move:
            return strprintf("move %zu", val);
        case command_type::pickup:
            return strprintf("pickup %zu", val);
        case command_type::charge_from_grid:
            return strprintf("charge_from_grid %zu", val);
        case command_type::charge_to_grid:
            return strprintf("charge_to_grid %zu", val);
        default:
            break;
        }
        return "";
    }
};
command parser(const string& command) {
    stringstream reader(command);
    string command_str;
    size_t value;
    reader >> command_str >> value;
    if (command_str == "stay") {
        return { command_type::stay, 0 };
    }
    else if (command_str == "move") {
        return { command_type::move, value };
    }
    else if (command_str == "pickup") {
        return { command_type::pickup, value };
    }
    else if (command_str == "charge_from_grid") {
        return { command_type::charge_from_grid, value };
    }
    else if (command_str == "charge_to_grid") {
        return { command_type::charge_to_grid, value };
    }
    return { invalid_command, (size_t)-1 };
}
int main() {
    setbuf(log_dest, nullptr);
    size_t N_solution = 1;
    cin >> N_solution;
    B prob(cin);
    std::shared_ptr<strategy<B>> str = nullptr;
    graph_summary gs(prob.graph, prob.grid);
    grid_info grid_i(prob.grid.N_grid);
    EV_info ev_i(prob.EV.N_EV);
    order_info order_i;
    string command_per_turn;
    vector<pair<double, double>> scores; scores.reserve(N_solution);
    for (size_t n = 0; n < N_solution; ++n) {
         //str.reset(new all_stay<B>(prob, gs));
        // str.reset(new random_walk<B>(prob, gs));
        // str.reset(new transport_only_0(prob, gs, n + 1));
        str.reset(new GreedyTransportStrategy(prob, gs, n + 1));
        str->initialize();

        set<int> picked;
        set<int> dropped;
        for (size_t t = 0; t < prob.T_max; ++t) {
            grid_i.load(cin);
            ev_i.load(cin);
            order_i.load(cin);
            str->command(grid_i, ev_i, order_i);
            command_per_turn = str->dequeue(ev_i);
            auto command_list = split_command(command_per_turn);
            cout << command_per_turn << flush;
            //cerr << command_per_turn << flush;

            {
                int left = 0;

                set<int> currentlyPicked;
                for (int i = 0; i < order_i.N_order; i++) {
                    if (!order_i.state[i]) {
                        left += 1;
                        continue;
                    }
                    currentlyPicked.insert(order_i.id[i]);
                }

                for (auto id : picked) {
                    if (currentlyPicked.count(id)) {
                        continue;
                    }
                    if (dropped.count(id)) {
                        cerr << "dup " << id << endl;
                    }
                    dropped.insert(id);
                }

                picked = currentlyPicked;

                if (t == prob.T_max - 1) {
                    cerr << left << "/" << picked.size() << "/" << dropped.size() << endl;
                }
            }
        }
        grid_i.load(cin);
        ev_i.load(cin);
        order_i.load(cin);
        double S_trans, S_ele;
        cin >> S_trans >> S_ele;
        cerr << S_trans << " " << S_ele << endl;
        cerr << (S_trans - prob.S_ref_trans) << " * " << (S_ele - prob.S_ref_ele) << " = " << (S_trans - prob.S_ref_trans) * (S_ele - prob.S_ref_ele) << endl;

    }
    return 0;
}
