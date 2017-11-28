package beast.evolution.datatype;

import beast.evolution.datatype.DataType.Base;
import beast.core.Description;
import beast.core.Input;

@Description("Datatype for sequences with twenty states V, G, A, D, P, E, S, L, R, T, I, Q, N, C, Y, H, M, W, K, F")
public class TwentyState extends Base {
	public Input<String> symbolOrderInput = new Input<String>("symbolOrder", "Order in which the symbols appear and their values (deafult VGEADPLR)");
    int[][] x = {
            {0},  // V
            {1},  // G
            {2},  // E
            {3},  // A
            {4},  // 
            {5},
            {6},
            {7},
            {8},
            {9},
            {10},
            {11},
            {12},
            {13},
            {14},
            {15},
            {16},
            {17},
            {18},
            {19},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}, // -
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19} // ?
    };

    public TwentyState() {
        stateCount = 20;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "VGADPESLRTIQNCYHMWKF" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public void initAndValidate() {
    	if (symbolOrderInput.get() != null) {
    		String symbolOrder = symbolOrderInput.get().trim();
    		symbolOrder = symbolOrder.replaceAll("[ \t\n]", "");
    		if (symbolOrder.length() != stateCount) {
    			throw new RuntimeException("symbolOrder should contain " + stateCount + " characters.");
    		}
    		codeMap = symbolOrder + GAP_CHAR + MISSING_CHAR;
    	}
    	
    	super.initAndValidate();
    }
    
    @Override
    public String getTypeDescription() {
        return "eightstate";
    }

}
