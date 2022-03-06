/**
 * 
 */
package nufeb.examples.ibm.surface;

import java.util.Map;

import eu.quanticol.jsstl.core.formula.AndFormula;
import eu.quanticol.jsstl.core.formula.AtomicFormula;
import eu.quanticol.jsstl.core.formula.EventuallyFormula;
import eu.quanticol.jsstl.core.formula.GloballyFormula;
import eu.quanticol.jsstl.core.formula.NotFormula;
import eu.quanticol.jsstl.core.formula.ParametricExpression;
import eu.quanticol.jsstl.core.formula.ParametricInterval;
import eu.quanticol.jsstl.core.formula.ReferencedFormula;
import eu.quanticol.jsstl.core.formula.SignalExpression;
import eu.quanticol.jsstl.core.formula.SomewhereFormula;
import eu.quanticol.jsstl.core.formula.jSSTLScript;
import eu.quanticol.jsstl.core.formula.EverywhereFormula;

/**
 * @author Bowen Li
 *
 */
public class IBMSurfaceScript extends jSSTLScript {

	// surface 
	public static final double minVF_CONST_ = 0.01;
	public static final double minVF_CONST_1 = 0.0;
	public static final int VF_VAR_ = 0;
	
	public IBMSurfaceScript() {
		super( 
			new String[] {
				"VF",
				"gridX"
			}
		);	
		
		// surface
		addFormula( "surface-s1" ,
				new AtomicFormula( 
						new ParametricExpression( ) {
						
							public SignalExpression eval( Map<String, Double> parameters ) {
								
								return new SignalExpression() {						
											
									public double eval(double... variables) {
										return variables[getIndex(VF_VAR_)] - minVF_CONST_;
									}	
														
								};	
											
							}
						
						} , 
				false
				)			
				 ,
		null );
		
		addFormula( "surface-s2" ,
				new SomewhereFormula( 
					new ParametricInterval( 
						new ParametricExpression() {
						
							public SignalExpression eval( final Map<String,Double> parameters ) {
					
								return new SignalExpression() {
									
									public double eval( double ... variables ) {
										return parameters.get("minD");
									}
									
								};					
								
							}
							
						} , 
						new ParametricExpression() {
						
							public SignalExpression eval( final Map<String,Double> parameters ) {
					
								return new SignalExpression() {
									
									public double eval( double ... variables ) {
										return parameters.get("maxD");
									}
									
								};					
								
							}
							
						} 		
					)		
					 ,
					new AtomicFormula( 
							new ParametricExpression( ) {
							
								public SignalExpression eval( Map<String, Double> parameters ) {
									
									return new SignalExpression() {						
												
										public double eval(double... variables) {
											return minVF_CONST_ - variables[getIndex(VF_VAR_)];
										}	
															
									};	
												
								}
							
							} , 
					false
					)	
				)		
				 ,
		null );
		
		addFormula( "surface" ,
				new AndFormula(
					new ReferencedFormula( 
						this ,
						"surface-s1"
						),
					new ReferencedFormula( 
						this ,
						"surface-s2"
					)
				),
		null );
		
		addFormula( "surface-find" ,
				new EventuallyFormula( 
						new ParametricInterval( 
							new ParametricExpression() {
							
								public SignalExpression eval( final Map<String,Double> parameters ) {
						
									return new SignalExpression() {
										
										public double eval( double ... variables ) {
											return parameters.get("minT");
										}
										
									};					
									
								}
								
							} , 
							new ParametricExpression() {
							
								public SignalExpression eval( final Map<String,Double> parameters ) {
						
									return new SignalExpression() {
										
										public double eval( double ... variables ) {
											return parameters.get("maxT");
										}
										
									};					
									
								}
								
							} 		
						),
						new ReferencedFormula( 
								this ,
								"surface"
					   )
					),
		null );
		
		addFormula( "layer-l1" ,
				new EverywhereFormula( 
					new ParametricInterval( 
						new ParametricExpression() {
						
							public SignalExpression eval( final Map<String,Double> parameters ) {
					
								return new SignalExpression() {
									
									public double eval( double ... variables ) {
										return 0;
									}
									
								};					
								
							}
							
						} , 
						new ParametricExpression() {
						
							public SignalExpression eval( final Map<String,Double> parameters ) {
					
								return new SignalExpression() {
									
									public double eval( double ... variables ) {
										return parameters.get("minD") - 1;
									}
									
								};					
								
							}
							
						} 		
					)		
					 ,
					new AtomicFormula( 
							new ParametricExpression( ) {
							
								public SignalExpression eval( Map<String, Double> parameters ) {
									
									return new SignalExpression() {						
												
										public double eval(double... variables) {
											return variables[getIndex(VF_VAR_)] - minVF_CONST_1;
										}	
															
									};	
												
								}
							
							} , 
					false
					)	
				)		
				 ,
		null );
		
		addFormula( "layer" ,
				new AndFormula(
					new ReferencedFormula( 
						this ,
						"surface"
						),
					new ReferencedFormula( 
						this ,
						"layer-l1"
					)
				),
		null );
		
		addFormula( "layer-find" ,
				new EventuallyFormula( 
						new ParametricInterval( 
							new ParametricExpression() {
							
								public SignalExpression eval( final Map<String,Double> parameters ) {
						
									return new SignalExpression() {
										
										public double eval( double ... variables ) {
											return parameters.get("minT");
										}
										
									};					
									
								}
								
							} , 
							new ParametricExpression() {
							
								public SignalExpression eval( final Map<String,Double> parameters ) {
						
									return new SignalExpression() {
										
										public double eval( double ... variables ) {
											System.out.println("maxT " + parameters.get("maxT"));
											return parameters.get("maxT");
										}
										
									};					
									
								}
								
							} 		
						),
						new ReferencedFormula( 
								this ,
								"layer"
					   )
					),
		null );
	}
}
