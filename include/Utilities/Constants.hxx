#ifndef T2S_CONSTANTS_HXX
#define T2S_CONSTANTS_HXX

namespace Tree2Sexaquark {

namespace Cuts {
enum Option_t {  //
    kMinimum,
    kMaximum,
    kAbsoluteMax
};
}  // namespace Cuts

enum Species_t {
    kLambda,
    kKaonZero,
    kPionPair,
};

enum Channel_t {
    kChannelA,
    kChannelD,
    kChannelE,
    kChannelH,
};

}  // namespace Tree2Sexaquark

#endif  // T2S_CONSTANTS_HXX
