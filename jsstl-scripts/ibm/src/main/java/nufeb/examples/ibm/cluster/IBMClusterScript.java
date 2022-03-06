/**
 * 
 */
package nufeb.examples.ibm.cluster;

import eu.quanticol.jsstl.core.formula.*;
import java.util.Map;		

public class IBMClusterScript extends jSSTLScript {

		public static final int VF_VAR_ = 0;
		public static final int gridX_VAR_ = 0;
		
		public static final double PRE_minVF_CONST_ = 0.01;
		public static final double POST_minVF_CONST = 0.02;
		
		public IBMClusterScript() {
			super( 
				new String[] {
					"VF",
					"gridX"
				}
			);	
			addFormula( "spot" ,
				new SurroundFormula( 
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
										return (variables[getIndex(VF_VAR_)] - PRE_minVF_CONST_);
									}	
														
								};	
											
							}
						
						} , 
						false
					)		
					 ,
					new AtomicFormula( 
						new ParametricExpression( ) {
						
							public SignalExpression eval( Map<String, Double> parameters ) {
								
								return new SignalExpression() {						
											
									public double eval(double... variables) {
										return (POST_minVF_CONST - (variables[getIndex(VF_VAR_)]));
									}	
														
								};	
											
							}
						
						} , 
						false
					)		
				)		
				 ,
				null );
			addFormula( "cluster-formation" ,
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
					)	 ,
					new ReferencedFormula( 
							this ,
							"spot"
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
