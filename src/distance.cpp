#include "distance.h"
#include "scorematrix.h"

namespace Distance{

bool USE_ENCODED_GAP = true;
char GAP = ScoreMatrix::GAP;
ScoreMatrix::EncodedResidue EGAP = ScoreMatrix::encode(GAP);

}
