#include "SystemLocant.h"
#include <iostream>
#include <vector>
#include <QString>

static int passed = 0;
static int failed = 0;

#define EXPECT_EQ(actual, expected) do { \
    auto actual_val_ = (actual); \
    auto expected_val_ = (expected); \
    if (actual_val_ != expected_val_) { \
        std::cerr << "FAIL at " << __LINE__ << ": expected " << #expected << " (" << expected_val_.toStdString() << "), got " << #actual << " (" << actual_val_.toStdString() << ")\n"; \
        failed++; \
    } else { \
        passed++; \
    } \
} while (0)

int main() {
    // TEST(SystemLocantTest, RendersIupacWithSuperscript)
    {
        SystemLocant loc{3, 'a', 1};
        EXPECT_EQ(renderIupac(loc), QString::fromUtf8("3a\u00B9"));
    }

    // TEST(SystemLocantTest, RendersCasLetterSuffix)
    {
        SystemLocant loc{10, 'b', 0};
        EXPECT_EQ(renderCas(loc), QString("10b"));
    }

    // TEST(SystemLocantTest, PyreneBothStylesFromSameInternalLocants)
    {
        SystemLocant a{3, 'a', 1}, b{5, 'a', 1};
        EXPECT_EQ(renderIupac(a) + "," + renderIupac(b), QString::fromUtf8("3a\u00B9,5a\u00B9"));
        SystemLocant casA{10, 'b', 0}, casB{10, 'c', 0};
        EXPECT_EQ(renderCas(casA) + "," + renderCas(casB), QString("10b,10c"));
    }

    // Additional edge cases:
    {
        SystemLocant plain{5, '\0', 0};
        EXPECT_EQ(renderIupac(plain), QString("5"));
        EXPECT_EQ(renderCas(plain), QString("5"));

        SystemLocant multiDigitSuperscript{3, 'a', 12};
        EXPECT_EQ(renderIupac(multiDigitSuperscript), QString::fromUtf8("3a\u00B9\u00B2"));
    }

    std::cout << "SystemLocant test results: " << passed << " passed, " << failed << " failed.\n";
    return failed == 0 ? 0 : 1;
}
