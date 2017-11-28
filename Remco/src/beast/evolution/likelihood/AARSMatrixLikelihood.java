package beast.evolution.likelihood;



import java.util.List;
import java.util.Random;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.evolution.substitutionmodel.AARSMatrixModel;

public class AARSMatrixLikelihood extends Distribution {
	public Input<AARSMatrixModel> matrixModelInput = new Input<AARSMatrixModel>("model","matrix model that generates indices arrays", Validate.REQUIRED);

	AARSMatrixModel matrixModel;
	
	@Override
	public void initAndValidate() {
		matrixModel = matrixModelInput.get();
	}
	
	@Override
	public double calculateLogP() {
		return Math.log(matrixModel.getFit());
	}
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
