$( document ).ready(function() {

  $(document).on('shiny:value', function(event) {
    if (event.name === 'AdaptativeSelection_ui_1-PredictiveFeatures' && event.value) {
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
              $('#placeholder4').append(
                '<div class=\"form-group shiny-input-container\">' +
                  '<label class=\"control-label\" id=\"AdaptativeSelection_ui_1-bdc' + i + '-label\" for=\"AdaptativeSelection_ui_1-bdc' + i + '\">Promoter length (from TSS)</label>' +
                  '<input id=\"AdaptativeSelection_ui_1-bdc' + i + '\" type=\"number\" class=\"form-control\" value=\"2000\"/>' +
                  '<br/>' +
                '</div>'
                
              );
              
              pwmname = 2000;
                
              } else if (myfeatures[i] === "ProximalPromoter") {
              $('#placeholder4').append(
                '<div class=\"form-group shiny-input-container\">' +
                  '<label class=\"control-label\" id=\"AdaptativeSelection_ui_1-bdc' + i + '-label\" for=\"AdaptativeSelection_ui_1-bdc' + i + '\">Proximal promoter length (from TSS)</label>' +
                  '<input id=\"AdaptativeSelection_ui_1-bdc' + i + '\" type=\"number\" class=\"form-control\" value=\"500\"/>' +
                  '<br/>' +
                '</div>'
                );
                
                pwmname = 500;
                
              } else if (myfeatures[i] === "Downstream") {
                
                $('#placeholder4').append(
                '<div class=\"form-group shiny-input-container\">' +
                  '<label class=\"control-label\" id=\"AdaptativeSelection_ui_1-bdc' + i + '-label\" for=\"AdaptativeSelection_ui_1-bdc' + i + '\">Downstream region length (from TTS)</label>' +
                  '<input id=\"AdaptativeSelection_ui_1-bdc' + i + '\" type=\"number\" class=\"form-control\" value=\"1000\"/>' +
                  '<br/>' +
                '</div>'
                );
                
                pwmname = 1000;
                
              }  else if (myfeatures[i] === "X5UTR") {
              $('#placeholder4').append(
                
              '<div class=\"form-group shiny-input-container\">' +
                '<label class=\"control-label\" id=\"AdaptativeSelection_ui_1-bdc' + i + '-label\"'  +
                          'for=\"AdaptativeSelection_ui_1-bdc' + i + '\">' +
                            '<div>' +
                              "5\'UTR" + ' data' +
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

                } else if (myfeatures[i] === "X3UTR") {
                                $('#placeholder4').append(
              '<div class=\"form-group shiny-input-container\">' +
                '<label class=\"control-label\" id=\"AdaptativeSelection_ui_1-bdc' + i + '-label\"'  +
                          'for=\"AdaptativeSelection_ui_1-bdc' + i + '\">' +
                            '<div>' +
                              "3\'UTR" + ' data' +
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

                } else {
                                $('#placeholder4').append(
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
                  
                


            } else {
                          $('#placeholder3').append(
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
