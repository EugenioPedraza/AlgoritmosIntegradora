#include "res.h"

std::string Algorithm::manacher(std::string transmissionText){ 
    int N = transmissionText.size();
    if (N == 0)
        return "";
    N = 2*N + 1; // Position count
    std::vector<int> L(N, 0); // LPS Length Array
    L[0] = 0;
    L[1] = 1;
    int C = 1; // centerPosition
    int R = 2; // centerRightPosition
    int i = 0; // currentRightPosition
    int iMirror; // currentLeftPosition
    int maxLPSLength = 0;
    int maxLPSCenterPosition = 0;
    int start = -1;
    int end = -1;
    int diff = -1;

    for (i = 2; i < N; i++) {
        // get currentLeftPosition iMirror for currentRightPosition i
        iMirror = 2*C-i;
        L[i] = 0;
        diff = R - i;
        // If currentRightPosition i is within centerRightPosition R
        if (diff > 0)
            L[i] = std::min(L[iMirror], diff);

        // Attempt to expand palindrome centered at currentRightPosition i
        while (((i + L[i]) < N && (i - L[i]) > 0) && 
            (((i + L[i] + 1) % 2 == 0) || 
            (transmissionText[(i + L[i] + 1)/2] == transmissionText[(i - L[i] - 1)/2]))) {
            L[i]++;
        }

        if (L[i] > maxLPSLength) { // Track maxLPSLength
            maxLPSLength = L[i];
            maxLPSCenterPosition = i;
        }

        // If palindrome centered at currentRightPosition i 
        // expand beyond centerRightPosition R,
        // adjust centerPosition C based on expanded palindrome.
        if (i + L[i] > R) {
            C = i;
            R = i + L[i];
        }
    }
    start = (maxLPSCenterPosition - maxLPSLength)/2;
    end = start + maxLPSLength - 1;
    return transmissionText.substr(start, end-start+1);
}