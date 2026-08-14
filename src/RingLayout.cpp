#include "RingLayout.h"

#include <cmath>

namespace {

RoutingOption makeRoutingOption(const Node& src, const Node& tgt, bool hvFirst) {
    RoutingOption opt;
    if (hvFirst) {
        if (src.x != tgt.x) opt.segments.push_back({src.x, src.y, tgt.x, src.y});
        if (src.y != tgt.y) opt.segments.push_back({tgt.x, src.y, tgt.x, tgt.y});
    } else {
        if (src.y != tgt.y) opt.segments.push_back({src.x, src.y, src.x, tgt.y});
        if (src.x != tgt.x) opt.segments.push_back({src.x, tgt.y, tgt.x, tgt.y});
    }
    return opt;
}

}  // namespace

RingLayout buildRingLayout(
        const std::vector<Node>& nodes,
        const MILPSolveResult& layout) {
    RingLayout out;
    out.tour = layout.tour;
    out.edges.reserve(layout.tourEdges.size());

    for (const auto& te : layout.tourEdges) {
        const Node& src = nodes[te.from];
        const Node& tgt = nodes[te.to];

        Edge e;
        e.from = te.from;
        e.to = te.to;
        e.distance = std::abs(src.x - tgt.x) + std::abs(src.y - tgt.y);
        e.option1 = makeRoutingOption(src, tgt, true);
        e.option2 = makeRoutingOption(src, tgt, false);
        out.edges.push_back(std::move(e));
    }

    out.routing.reserve(layout.tourEdges.size());
    for (size_t i = 0; i < layout.tourEdges.size(); ++i) {
        const auto& te = layout.tourEdges[i];
        const Edge& stored = out.edges[i];
        int optIdx = te.routing == 0 ? 1 : 2;
        const RoutingOption* opt = optIdx == 1 ? &stored.option1 : &stored.option2;
        out.routing.push_back({&stored, optIdx, opt});
    }

    return out;
}

std::vector<Segment> pathVerticesToSegments(
        const std::vector<std::pair<double, double>>& vertices) {
    std::vector<Segment> segs;
    for (size_t i = 0; i + 1 < vertices.size(); ++i) {
        segs.push_back({
            vertices[i].first, vertices[i].second,
            vertices[i + 1].first, vertices[i + 1].second});
    }
    return segs;
}
