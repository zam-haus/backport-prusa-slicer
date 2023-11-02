#include "libslic3r/Point.hpp"
#include <catch2/catch_all.hpp>
#include <libslic3r/SupportSpotsGenerator.hpp>

using namespace Slic3r;
using namespace SupportSpotsGenerator;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


TEST_CASE("Numerical integral calculation compared with exact solution.", "[SupportSpotsGenerator]") {
    const float width = 10;
    const float height = 20;
    const Polygon polygon = {
        scaled(Vec2f{-width / 2, -height / 2}),
        scaled(Vec2f{width / 2, -height / 2}),
        scaled(Vec2f{width / 2, height / 2}),
        scaled(Vec2f{-width / 2, height / 2})
    };

    const Integrals integrals{{polygon}};
    CHECK_THAT(integrals.area, WithinRel(width * height));
    CHECK_THAT(integrals.x_i.x(), WithinRel(0.));
    CHECK_THAT(integrals.x_i.y(), WithinRel(0.));
    CHECK_THAT(integrals.x_i_squared.x(), WithinRel(std::pow(width, 3) * height / 12, EPSILON));
    CHECK_THAT(integrals.x_i_squared.y(), WithinRel(width * std::pow(height, 3) / 12, EPSILON));
}

TEST_CASE("Moment values and ratio check.", "[SupportSpotsGenerator]") {
    const float width = 40;
    const float height = 2;

    // Moments are calculated at centroid.
    // Polygon centroid must not be (0, 0).
    const Polygon polygon = {
        scaled(Vec2f{0, 0}),
        scaled(Vec2f{width, 0}),
        scaled(Vec2f{width, height}),
        scaled(Vec2f{0, height})
    };

    const Integrals integrals{{polygon}};

    const Vec2f x_axis{1, 0};
    const float x_axis_moment = compute_second_moment(integrals, x_axis);

    const Vec2f y_axis{0, 1};
    const float y_axis_moment = compute_second_moment(integrals, y_axis);

    const float moment_ratio = std::pow(width / height, 2);

    // Ensure the object transaltion has no effect.
    CHECK_THAT(x_axis_moment, WithinRel(width * std::pow(height, 3) / 12, EPSILON));
    CHECK_THAT(y_axis_moment, WithinRel(std::pow(width, 3) * height / 12, EPSILON));
    // If the object is "wide" the y axis moments should be large compared to x axis moment.
    CHECK_THAT(y_axis_moment / x_axis_moment, WithinRel(moment_ratio));
}

TEST_CASE("Moments calculation for rotated axis.", "[SupportSpotsGenerator]") {

    Polygon polygon = {
        scaled(Vec2f{6.362284076172198, 138.9674202217155}),
        scaled(Vec2f{97.48779843751677, 106.08136606617076}),
        scaled(Vec2f{135.75221821532384, 66.84428834668765}),
        scaled(Vec2f{191.5308049852741, 45.77905628725614}),
        scaled(Vec2f{182.7525148049201, 74.01799041087513}),
        scaled(Vec2f{296.83210979283473, 196.80022572637228}),
        scaled(Vec2f{215.16434429179148, 187.45715418834143}),
        scaled(Vec2f{64.64574271229334, 284.293883209721}),
        scaled(Vec2f{110.76507036894843, 174.35633141113783}),
        scaled(Vec2f{77.56229640885199, 189.33057746591336})
    };

    Integrals integrals{{polygon}};

    // Meassured counterclockwise from (1, 0)
    const float angle = 1.432f;
    Vec2f axis{std::cos(angle), std::sin(angle)};

    float moment_calculated_then_rotated = compute_second_moment(
        integrals,
        axis
    );

    // We want to rotate the object clockwise by angle to align the axis with (1, 0)
    // Method .rotate is counterclockwise for positive angle
    polygon.rotate(-angle);

    Integrals integrals_rotated{{polygon}};
    float moment_rotated_polygon = compute_second_moment(
        integrals_rotated,
        Vec2f{1, 0}
    );

    // Up to 0.1% accuracy
    CHECK_THAT(moment_calculated_then_rotated, Catch::Matchers::WithinRel(moment_rotated_polygon, 0.001f));
}

struct ObjectPartFixture {
    const Polyline polyline{
        Point{scaled(Vec2f{0, 0})},
        Point{scaled(Vec2f{1, 0})},
    };
    const float width = 0.1f;
    bool connected_to_bed = true;
    coordf_t print_head_z = 0.2;
    coordf_t layer_height = 0.2;
    ExtrusionAttributes attributes;
    ExtrusionEntityCollection collection;
    std::vector<const ExtrusionEntityCollection*> extrusions{};
    Polygon expected_polygon{
        Point{scaled(Vec2f{0, -width / 2})},
        Point{scaled(Vec2f{1, -width / 2})},
        Point{scaled(Vec2f{1, width / 2})},
        Point{scaled(Vec2f{0, width / 2})}
    };

    ObjectPartFixture() {
        attributes.width = width;
        const ExtrusionPath path{polyline, attributes};
        collection.append(path);
        extrusions.push_back(&collection);
    }
};

TEST_CASE_METHOD(ObjectPartFixture, "Constructing ObjectPart using extrusion collections", "[SupportSpotsGenerator]") {
    ObjectPart part{
        extrusions,
        connected_to_bed,
        print_head_z,
        layer_height,
        std::nullopt
    };

    Integrals expected{{expected_polygon}};

    CHECK(part.connected_to_bed == true);
    Vec3f volume_centroid{part.volume_centroid_accumulator / part.volume};
    CHECK_THAT(volume_centroid.x(), WithinAbs(0.5, EPSILON));
    CHECK_THAT(volume_centroid.y(), WithinAbs(0., EPSILON));
    CHECK_THAT(volume_centroid.z(), WithinAbs(layer_height / 2, EPSILON));
    CHECK_THAT(part.sticking_area, WithinRel(expected.area));
    CHECK_THAT(part.sticking_centroid_accumulator.x(), WithinRel(expected.x_i.x()));
    CHECK_THAT(part.sticking_centroid_accumulator.y(), WithinRel(expected.x_i.y()));
    CHECK_THAT(part.sticking_second_moment_of_area_accumulator.x(), WithinRel(expected.x_i_squared.x()));
    CHECK_THAT(part.sticking_second_moment_of_area_accumulator.y(), WithinRel(expected.x_i_squared.y()));
    CHECK_THAT(part.sticking_second_moment_of_area_covariance_accumulator, WithinAbs(expected.xy, 1e-6));
    CHECK_THAT(part.volume, WithinRel(layer_height * width, EPSILON));
}

TEST_CASE_METHOD(ObjectPartFixture, "Constructing ObjectPart with brim", "[SupportSpotsGenerator]") {
    float brim_width = 1;
    Polygons brim = get_brim(ExPolygon{expected_polygon}, BrimType::btOuterOnly, brim_width);

    ObjectPart part{
        extrusions,
        connected_to_bed,
        print_head_z,
        layer_height,
        brim
    };

    CHECK_THAT(part.sticking_area, WithinRel((1 + 2*brim_width) * (width + 2*brim_width)));
}

