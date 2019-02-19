package beast.evolution.datatype;

import beast.evolution.datatype.DataType.Base;
import beast.core.Description;
import beast.core.Input;

@Description("Datatype for sequences with eight states V, G, E, A, D, P, L, R")
public class EightState extends Base {
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
            {0, 1, 2, 3, 4, 5, 6, 7}, // -
            {0, 1, 2, 3, 4, 5, 6, 7} // ?
    };

    public EightState() {
        stateCount = 8;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "VGEADPLR" + GAP_CHAR + MISSING_CHAR;
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
