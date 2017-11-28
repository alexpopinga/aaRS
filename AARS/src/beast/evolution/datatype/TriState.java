package beast.evolution.datatype;

import beast.evolution.datatype.DataType.Base;
import beast.core.Description;

@Description("Datatype for sequences with three states E,G,V")
public class TriState extends Base {
    int[][] x = {
            {0},  // V
            {1},  // E
            {2},  // G
            {0, 1, 2}, // -
            {0, 1, 2} // ?
    };

    public TriState() {
        stateCount = 3;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "VEG" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "tristate";
    }

}
