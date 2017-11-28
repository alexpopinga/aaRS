package beast.evolution.datatype;

import beast.evolution.datatype.DataType.Base;
import beast.core.Description;
import beast.core.Input;

@Description("Datatype for sequences with twenty states V, G, A, D, P, E, S, L, R, T, I, N, C, Y, H, M, W, K, F")
public class NineTeenState extends Base {
	public Input<String> symbolOrderInput = new Input<String>("symbolOrder", "Order in which the symbols appear and their values (deafult VGADPESLRTINCYHMWKF)");
    int[][] x = {
            {0},  
            {1},  
            {2},  
            {3},  
            {4},  
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
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}, // -
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18} // ?
    };

    public NineTeenState() {
        stateCount = 19;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "VGADPESLRTINCYHMWKF" + GAP_CHAR + MISSING_CHAR;
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
        return "nineteenstate";
    }

}
