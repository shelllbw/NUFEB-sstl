/**
 * 
 */
package nufeb.examples.ibm.shear;

import java.util.Map;

import eu.quanticol.jsstl.core.formula.AndFormula;
import eu.quanticol.jsstl.core.formula.AtomicFormula;
import eu.quanticol.jsstl.core.formula.EventuallyFormula;
import eu.quanticol.jsstl.core.formula.GloballyFormula;
import eu.quanticol.jsstl.core.formula.ParametricExpression;
import eu.quanticol.jsstl.core.formula.ParametricInterval;
import eu.quanticol.jsstl.core.formula.ReferencedFormula;
import eu.quanticol.jsstl.core.formula.SignalExpression;
import eu.quanticol.jsstl.core.formula.jSSTLScript;

/**
 * @author Bowen Li
 *
 */
public class IBMShearScript extends jSSTLScript {

	public static final int VF_VAR_ = 0;
	public static final int gridX_VAR_ = 0;
	
	public IBMShearScript() {
		super( 
			new String[] {
				"VF",
				"gridX"
			}
		);	
		
		// streamer formation F[O, maxT]G[0, 10h] VF > 0
		addFormula( "streamer-global" ,
				new GloballyFormula( 
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
											return 86400 ;
										}
										
									};					
									
								}
								
							} 		
						),
						new AtomicFormula( 
								new ParametricExpression( ) {
								
									public SignalExpression eval( Map<String, Double> parameters ) {
										
										return new SignalExpression() {						
													
											public double eval(double... variables) {
												return variables[getIndex(VF_VAR_)];
											}									
										};													
									}							
								} , 
						true
						)	
					),
		null );
		
		addFormula( "streamer-formation" ,
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
											//System.out.println("maxT " + parameters.get("maxT"));
											return parameters.get("maxT");
										}
										
									};					
									
								}
								
							} 		
						),
						new ReferencedFormula( 
								this ,
								"streamer-global"
						)
					),
		null );
		
		// detachment
		
		addFormula( "detach-find" ,
				new EventuallyFormula( 
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
											//System.out.println("maxT " + parameters.get("maxT"));
											return 7200;
										}
										
									};					
									
								}
								
							} 		
						),
						new AtomicFormula( 
								new ParametricExpression( ) {
								
									public SignalExpression eval( Map<String, Double> parameters ) {
										
										return new SignalExpression() {						
													
											public double eval(double... variables) {
												return 0.005 - variables[getIndex(VF_VAR_)];
											}									
										};													
									}							
								} , 
						true
						)	
					),
		null );
		
		addFormula( "detach-and" ,
			new AndFormula(
					new AtomicFormula( 
							new ParametricExpression( ) {
							
								public SignalExpression eval( Map<String, Double> parameters ) {
									
									return new SignalExpression() {						
												
										public double eval(double... variables) {
											return variables[getIndex(VF_VAR_)];
										}									
									};													
								}							
							} , 
					true
					),
					new ReferencedFormula( 
							this ,
							"detach-find"
					)
			),
		null );
		
		addFormula( "detachment" ,
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
								"detach-and"
						)	
					),
		null );
		
		// grid bound GridX > fromX && GridX < toX
		addFormula( "grid-bound" ,
				new GloballyFormula( 
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
											return parameters.get("maxT");
										}
										
									};					
									
								}
								
							} 		
						),
						new AndFormula(
								new AtomicFormula( 
										new ParametricExpression( ) {
										
											public SignalExpression eval(final Map<String, Double> parameters ) {
												
												return new SignalExpression() {						
															
													public double eval(double... variables) {
														return variables[getIndex(gridX_VAR_)] - parameters.get("fromX");
													}	
																		
												};	
															
											}
										
										} , 
								false
								),
								new AtomicFormula( 
										new ParametricExpression( ) {
										
											public SignalExpression eval(final Map<String, Double> parameters ) {
												
												return new SignalExpression() {						
															
													public double eval(double... variables) {
														return  parameters.get("toX") - variables[getIndex(gridX_VAR_)];
													}	
																		
												};	
															
											}
										
										} , 
								false
								)	
						)
					),
		null );
	}
}
