#include <math.h>
#include <vector>
#include <chrono>

using namespace std;

float fastSqrt(float x) {
    union {
        int32_t i;
        float x;
    } u;

    u.x = x;
    u.i = (1 << 29) + (u.i >> 1) - (1 << 22);
    return u.x;
}

/// Helper structure to for measuring the time taken by execution of different functions
template <typename TimeT = chrono::milliseconds, typename ClockT = chrono::steady_clock>
struct measure {
    template <typename F, typename... Args>
    static TimeT execution(F func, Args&&... args) {
        auto start = ClockT::now();

        // The function invocation we want to measure
        func(forward<Args>(args)...);

        return chrono::duration_cast<TimeT>(ClockT::now() - start);
    }
};

namespace Warmup {

void skipTraverse() {
    constexpr int numValues = 1000'000;
    constexpr int repeatCount = 5'000;

    printf("Generating input values...\n");
    fflush(stdout);

    // Generate a lot of random numbers
    vector<int> values;
    values.reserve(numValues);
    for (int i = 0; i < numValues; i++)
        values.push_back(rand());

    printf("Performing the skipTraverse tests...\n\n");
    fflush(stdout);

    // Measure time for incrementing
    for (int step = 1; step < 130; step++) {
        auto durMs = measure<>::execution([&] {
            for (int r = 0; r < repeatCount; r++) {
                for (int i = 0; i < values.size(); i += step)
                    values[i]++;
            }
        });

        printf("step=%4d, NumOper=%7d, stride=%5d, time=%5d ms\n", step,
                int(values.size() / step), int(step * sizeof(int)), int(durMs.count()));
        fflush(stdout);
    }

    printf("\n\n");
    fflush(stdout);

    // Square each element - 1 multiplication
    for (int step = 4; step < 17; step++) {
        auto durMs = measure<>::execution([&] {
            for (int r = 0; r < repeatCount; r++) {
                for (int i = 0; i < values.size()-1; i += step)
                    values[i] = values[i]*values[i];
            }
        });

        printf("Mul1 - step=%4d, NumOper=%7d, stride=%5d, time=%5d ms\n", step,
                int(values.size() / step), int(step * sizeof(int)), int(durMs.count()));
        fflush(stdout);
    }

    printf("\n\n");
    fflush(stdout);

    // Square each element - 2 multiplications
    for (int step = 4; step < 17; step++) {
        auto durMs = measure<>::execution([&] {
            for (int r = 0; r < repeatCount; r++) {
                for (int i = 0; i < values.size()-1; i += step)
                    values[i] = values[i]*values[i]*values[i];
            }
        });

        printf("Mul2 - step=%4d, NumOper=%7d, stride=%5d, time=%5d ms\n", step,
                int(values.size() / step), int(step * sizeof(int)), int(durMs.count()));
        fflush(stdout);
    }

    printf("\n\n");
    fflush(stdout);

    // Square each element - 3 multiplications
    for (int step = 4; step < 17; step++) {
        auto durMs = measure<>::execution([&] {
            for (int r = 0; r < repeatCount; r++) {
                for (int i = 0; i < values.size()-1; i += step)
                    values[i] = values[i]*values[i]*values[i]*values[i];
            }
        });

        printf("Mul3 - step=%4d, NumOper=%7d, stride=%5d, time=%5d ms\n", step,
                int(values.size() / step), int(step * sizeof(int)), int(durMs.count()));
        fflush(stdout);
    }

    printf("\n\n");
    fflush(stdout);

    // Square each element - 3 multiplications
    for (int step = 4; step < 17; step++) {
        auto durMs = measure<>::execution([&] {
            for (int r = 0; r < repeatCount; r++) {
                for (int i = 0; i < values.size()-1; i += step)
                    values[i] = values[i]*values[i]*values[i]*values[i]*values[i];
            }
        });

        printf("Mul4 - step=%4d, NumOper=%7d, stride=%5d, time=%5d ms\n", step,
                int(values.size() / step), int(step * sizeof(int)), int(durMs.count()));
        fflush(stdout);
    }

    printf("\n\n");
    fflush(stdout);

    // Now add a suqare root
    for (int step = 4; step < 17; step++) {
        auto durMs = measure<>::execution([&] {
            for (int r = 0; r < repeatCount; r++) {
                for (int i = 0; i < values.size(); i += step)
                    values[i] = fastSqrt(values[i]);
            }
        });

        printf("FastSqrt - step=%4d, NumOper=%7d, stride=%5d, time=%5d ms\n", step,
                int(values.size() / step), int(step * sizeof(int)), int(durMs.count()));
        fflush(stdout);
    }
}

} // namespace Warmup

namespace SqrtTest {

void testSpeed() {
    constexpr long repeatCount = 1'000'000'000;

    float sum = 0.0;

    // Measure the time needed to perform just the plumbing code
    auto durMsBase = measure<>::execution([&] {
        for (long r = 0; r < repeatCount; r++) {
            sum += float(r);
        }
    });

    // Measure the time needed to perform a lot of sqrtf (plus some plumbing code)
    auto durMsSqrt = measure<>::execution([&] {
        for (long r = 0; r < repeatCount; r++) {
            sum += sqrtf(float(r));
        }
    });

    // Measure the time needed to perform a lot of our sqrt (plus some plumbing code)
    auto durMsFastSqrt = measure<>::execution([&] {
        for (long r = 0; r < repeatCount; r++) {
            sum += fastSqrt(float(r));
        }
    });

    if (sum < 1.0)
        printf(".");

    auto scale = 1'000'000 / double(repeatCount);
    auto diff1 = double(durMsSqrt.count() - durMsBase.count()) * scale;
    auto diff2 = double(durMsFastSqrt.count() - durMsBase.count()) * scale;
    printf("Cost of sqrtf: approx %.2f ns\n", diff1);
    printf("Cost of FastSqrt: approx %.2f ns\n", diff2);
    printf("   - sqrtf + plumbing: %.2f ns\n", durMsSqrt.count() * scale);
    printf("   - FastSqrt + plumbing: %.2f ns\n", durMsFastSqrt.count() * scale);
    printf("   - plumbing only: %.2f ns\n", durMsBase.count() * scale);
    fflush(stdout);
}

} // namespace SqrtTest

namespace ConcatProblem {
struct Vec3 {
    double x;
    double y;
    double z;
};

struct Vec3Compact {
    float x;
    float y;

    Vec3Compact() {}
    Vec3Compact(float xx, float yy)
        : x{xx}
        , y{yy} {}
    Vec3Compact(Vec3 v)
        : x{(float)v.x}
        , y{(float)v.y} {}

    explicit operator Vec3() const { return Vec3{x, y, 1.0 - fastSqrt(x * x + y * y)}; }
};

auto compact(const vector<Vec3>& src) {
    vector<Vec3Compact> res{src.begin(), src.end()};
    return res;
}

double genRand(double max = 1.0) { return max * (rand() % 10'000) / 10'000.0; }

vector<Vec3> generateRnd(int count) {
    vector<Vec3> res;
    res.resize(count);
    for (int i = 0; i < count; i++) {
        Vec3 cur;
        cur.x = genRand(1.0);
        cur.y = genRand(1.0 - cur.x);
        cur.z = 1.0 - sqrt(cur.x * cur.x + cur.y * cur.y);
        res[i] = cur;
    }
    return res;
}

void compare(const vector<Vec3>& lhs, const vector<Vec3>& rhs) {
    if (lhs.size() != rhs.size())
        throw runtime_error("Sizes don't match");
    for (int i = 0; i < lhs.size(); i++) {
        Vec3 val1 = lhs[i];
        Vec3 val2 = rhs[i];
        double dx = val2.x - val1.x;
        double dy = val2.y - val1.y;
        double dz = val2.z - val1.z;
        constexpr double errLimit = 0.1;
        if (fabs(dx) > errLimit || fabs(dy) > errLimit || fabs(dz) > errLimit) {
            printf("%d: (%.2f, %.2f, %.2f) != (%.2f, %.2f, %.2f)\n", i, val1.x, val1.y, val1.z,
                    val2.x, val2.y, val2.z);
            throw runtime_error("values don't match");
        }
    }
}

Vec3 negateXY(Vec3 v) { return Vec3{-v.x, -v.y, v.z}; }
Vec3 rotateZ30(Vec3 v) {
    constexpr float cos_theta = 0.866025f;
    constexpr float sin_theta = 0.5f;
    return Vec3{v.x * cos_theta - v.y * sin_theta, v.x * sin_theta + v.y * cos_theta, v.z};
}

Vec3 negateXY_tr(Vec3Compact v) { return (Vec3)Vec3Compact{-v.x, -v.y}; }
Vec3 rotateZ30_tr(Vec3Compact v) {
    constexpr float cos_theta = 0.866025f;
    constexpr float sin_theta = 0.5f;
    return (Vec3)Vec3Compact{v.x * cos_theta - v.y * sin_theta, v.x * sin_theta + v.y * cos_theta};
}

void concat_memcpy(
        const vector<Vec3>& s1, const vector<Vec3>& s2, vector<Vec3>& out, int repeatCount) {
    int sz1 = int(s1.size());
    int sz2 = int(s2.size());
    for (int r = 0; r < repeatCount; r++) {
        memcpy(&out[0], &s1[0], sz1 * sizeof(Vec3));
        memcpy(&out[sz1], &s2[0], sz2 * sizeof(Vec3));
    }
}

void concat_manual(
        const vector<Vec3>& s1, const vector<Vec3>& s2, vector<Vec3>& out, int repeatCount) {
    int sz1 = int(s1.size());
    int sz2 = int(s2.size());
    for (int r = 0; r < repeatCount; r++) {
        int idxDst = 0;
        for (int i = 0; i < sz1; i++, idxDst++)
            out[idxDst] = s1[i];
        for (int i = 0; i < sz2; i++, idxDst++)
            out[idxDst] = s2[i];
    }
}

void concat_andNegate(
        const vector<Vec3>& s1, const vector<Vec3>& s2, vector<Vec3>& out, int repeatCount) {
    int sz1 = int(s1.size());
    int sz2 = int(s2.size());
    for (int r = 0; r < repeatCount; r++) {
        int idxDst = 0;
        for (int i = 0; i < sz1; i++, idxDst++)
            out[idxDst] = negateXY(s1[i]);
        for (int i = 0; i < sz2; i++, idxDst++)
            out[idxDst] = negateXY(s2[i]);
    }
}

void concat_andRotate(
        const vector<Vec3>& s1, const vector<Vec3>& s2, vector<Vec3>& out, int repeatCount) {
    int sz1 = int(s1.size());
    int sz2 = int(s2.size());
    for (int r = 0; r < repeatCount; r++) {
        int idxDst = 0;
        for (int i = 0; i < sz1; i++, idxDst++)
            out[idxDst] = rotateZ30(s1[i]);
        for (int i = 0; i < sz2; i++, idxDst++)
            out[idxDst] = rotateZ30(s2[i]);
    }
}

void concat_compact(const vector<Vec3Compact>& s1, const vector<Vec3Compact>& s2, vector<Vec3>& out,
        int repeatCount) {
    int sz1 = int(s1.size());
    int sz2 = int(s2.size());
    for (int r = 0; r < repeatCount; r++) {
        int idxDst = 0;
        for (int i = 0; i < sz1; i++, idxDst++)
            out[idxDst] = (Vec3)s1[i];
        for (int i = 0; i < sz2; i++, idxDst++)
            out[idxDst] = (Vec3)s2[i];
    }
}

void concat_compact_andNegate(const vector<Vec3Compact>& s1, const vector<Vec3Compact>& s2,
        vector<Vec3>& out, int repeatCount) {
    int sz1 = int(s1.size());
    int sz2 = int(s2.size());
    for (int r = 0; r < repeatCount; r++) {
        int idxDst = 0;
        for (int i = 0; i < sz1; i++, idxDst++)
            out[idxDst] = negateXY_tr(s1[i]);
        for (int i = 0; i < sz2; i++, idxDst++)
            out[idxDst] = negateXY_tr(s2[i]);
    }
}

void concat_compact_andRotate(const vector<Vec3Compact>& s1, const vector<Vec3Compact>& s2,
        vector<Vec3>& out, int repeatCount) {
    int sz1 = int(s1.size());
    int sz2 = int(s2.size());
    for (int r = 0; r < repeatCount; r++) {
        int idxDst = 0;
        for (int i = 0; i < sz1; i++, idxDst++)
            out[idxDst] = rotateZ30_tr(s1[i]);
        for (int i = 0; i < sz2; i++, idxDst++)
            out[idxDst] = rotateZ30_tr(s2[i]);
    }
}

template <int skip>
void concat_compact_andRotateSkip(const vector<Vec3Compact>& s1, const vector<Vec3Compact>& s2,
        vector<Vec3>& out, int repeatCount) {
    int sz1 = int(s1.size());
    int sz2 = int(s2.size());
    for (int r = 0; r < repeatCount; r++) {
        int idxDst = 0;
        for (int i = 0; i < sz1; i++, idxDst++)
            if (i % skip == 0)
                out[idxDst] = rotateZ30_tr(s1[i]);
            else
                out[idxDst] = (Vec3)s1[i];
        for (int i = 0; i < sz2; i++, idxDst++)
            if (i % skip == 0)
                out[idxDst] = rotateZ30_tr(s2[i]);
            else
                out[idxDst] = (Vec3)s2[i];
    }
}

void doTest() {
    constexpr int numValues = 1000'000;
    constexpr int repeatCount = 100;

    // Original source vectors
    auto src1 = generateRnd(numValues);
    auto src2 = generateRnd(numValues);

    // The vector where we put the concatenation result
    vector<Vec3> res1, res2, res3, res1c, res2c, res3c, res4c;
    res1.resize(numValues * 2);
    res2.resize(numValues * 2);
    res3.resize(numValues * 2);
    res1c.resize(numValues * 2);
    res2c.resize(numValues * 2);
    res3c.resize(numValues * 2);
    res4c.resize(numValues * 2);

    // Simple concatenation - memcpy
    auto durMs = measure<>::execution(concat_memcpy, src1, src2, res1, repeatCount);
    printf("memcpy concat\ttime: %d ms\n", (int)durMs.count());
    fflush(stdout);

    // Simple concatenation -- second time, manual
    durMs = measure<>::execution(concat_manual, src1, src2, res1, repeatCount);
    printf("Simple concat\ttime: %d ms\n", (int)durMs.count());
    fflush(stdout);

    // Concatenation & negate
    durMs = measure<>::execution(concat_andNegate, src1, src2, res2, repeatCount);
    printf("Concat & negate\ttime: %d ms\n", (int)durMs.count());
    fflush(stdout);

    // Concatenation & rotate
    durMs = measure<>::execution(concat_andRotate, src1, src2, res3, repeatCount);
    printf("Concat & rotate\ttime: %d ms\n", (int)durMs.count());
    fflush(stdout);

    // Compact versions of those vectors
    auto src1c = compact(src1);
    auto src2c = compact(src2);

    // Compacted: simple concatenation
    durMs = measure<>::execution(concat_compact, src1c, src2c, res1c, repeatCount);
    printf("Compact: simple concat\ttime: %d ms\n", (int)durMs.count());
    fflush(stdout);

    // Compacted: concatenation & negation
    durMs = measure<>::execution(concat_compact_andNegate, src1c, src2c, res2c, repeatCount);
    printf("Compact: concat & negate\ttime: %d ms\n", (int)durMs.count());
    fflush(stdout);

    // Compacted: concatenation & rotation
    durMs = measure<>::execution(concat_compact_andRotate, src1c, src2c, res3c, repeatCount);
    printf("Compact: concat & rotate\ttime: %d ms\n", (int)durMs.count());
    fflush(stdout);

    // Compacted: concatenation & rotation-skip
    durMs = measure<>::execution(concat_compact_andRotateSkip<8>, src1c, src2c, res4c, repeatCount);
    printf("Compact: concat & rotate-skip\ttime: %d ms\n", (int)durMs.count());
    fflush(stdout);

    // Ensure correctness
    compare(res1, res1c);
    compare(res2, res2c);
    compare(res3, res3c);

    printf("Done.\n");
    fflush(stdout);
}

} // namespace ConcatProblem

int main(int argc, char** argv) {
    srand(1);

    Warmup::skipTraverse();
    ConcatProblem::doTest();

    // SqrtTest::testSpeed();

    return 0;
}
