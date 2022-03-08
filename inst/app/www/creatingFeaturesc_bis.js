$( document ).ready(function() {
  
  $(document).on('shiny:value', function(event) {
    if (event.name === 'AdaptativeSelection_ui_1-PredictiveFeaturesc' && event.value) {
      prevText = event.value;
      var pwmname = "";
      var myfeatures = event.value.split("\"The predictive features are: ");
      myfeatures = myfeatures[1];
      myfeatures = myfeatures.split(", ");
      myfeatures[myfeatures.length-1] = myfeatures[myfeatures.length-1].split("\"")[0];
      console.log("test");
      console.log(myfeatures);
      Shiny.unbindAll();
      if (myfeatures[0] !== "") {
        for (let i = 0; i < myfeatures.length; i++){
          if (["Promoter", "ProximalPromoter", "X5UTR", "Intron", "CDS", "X3UTR", "Downstream"].includes(myfeatures[i])) {
            if (myfeatures[i] === "Promoter"){

              pwmname = 2000;
              
            } else if (myfeatures[i] === "ProximalPromoter") {

              pwmname = 500;
              
            } else if (myfeatures[i] === "Downstream") {
              

              pwmname = 1000;
              
            }  else if (myfeatures[i] === "X5UTR") {

              pwmname = "";
              
            } else if (myfeatures[i] === "X3UTR") {

              pwmname = "";
              
            } else {

              pwmname = "";
              
            }
            
            
            
            
          } else {
            $('#placeholder5').append(
              '<div class=\"form-group shiny-input-container\">' +
                '<label class=\"control-label\" id=\"AdaptativeSelection_ui_1-bdc' + i + '-label\"'  +
                'for=\"AdaptativeSelection_ui_1-bdc' + i + '\">' +
                '<div>' +
                myfeatures[i] + ' data' +
                '</div>' +
                '</label>' +
                '<div class=\"input-group\">' +
                '<label class=\"input-group-btn input-group-prepend\">' +
                '<span class=\"btn btn-default btn-file\">' +
                'Browse...' +
                '<input id=\"AdaptativeSelection_ui_1-bdc' + i +  '\" name=\"AdaptativeSelection_ui_1-bdc' + i + '\" type=\"file\"' +
                'style=\"position: absolute !important; top: -99999px !important; left: -99999px !important;\"' +
                ' multiple=\"multiple\" accept=\".bed, .gtf\"/>' +
                '</span>' +
                '</label>' +
                '<input type=\"text\" class=\"form-control\" placeholder=\"No file(s) selected\" readonly=\"readonly\"/>' +
                '</div>' +
                '</div>'
            );
            
            pwmname = "";
            
          }
          
          $('#AdaptativeSelection_ui_1-bdc'+ i).val(pwmname);
          $('#AdaptativeSelection_ui_1-bdc'+ i).type = 'shiny:input';
        }
      }
      Shiny.bindAll();
    }
  });
});
