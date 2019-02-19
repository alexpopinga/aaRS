package beast.evolution.datatype;

import beast.evolution.datatype.DataType.Base;
import beast.core.Description;
import beast.core.Input;

@Description("Datatype for sequences with four states V,G,A,E")
public class FourState extends Base {
	public Input<String> symbolOrderInput = new Input<String>("symbolOrder", "Order in which the symbols appear and their values (deafult VGAE)");
    int[][] x = {
            {0},  // V
            {1},  // G
            {2},  // A
            {3},  // E
            {0, 1, 2, 3}, // -
            {0, 1, 2, 3} // ?
    };

    public FourState() {
        stateCount = 4;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "VGAE" + GAP_CHAR + MISSING_CHAR;
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
        return "fourstate";
    }

}
