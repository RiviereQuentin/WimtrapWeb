$( document ).ready(function() {

  $(document).on('shiny:inputchanged', function(event) {
    if (event.name === 'AdaptativeSelection_ui_1-filesData' && event.value) {
      prevText = event.value;
       Shiny.unbindAll();
      for (let i = 0; i < event.value.length; i++){
        $('#placeholder2').append(
        '<div class=\"form-group shiny-input-container\">' +
          '<label for=\"AdaptativeSelection_ui_1-bed' + i + '\">' + event.value[i].name +  '</label>' +
          '<input id=\"AdaptativeSelection_ui_1-bed' + i + '\" type=\"text\"' +
                  'class=\"form-control\" value=\"\">' + 
        '</div>'
        );
        var pwmname = "";
        $('#AdaptativeSelection_ui_1-bed'+ i).val(pwmname);
        $('#AdaptativeSelection_ui_1-bed'+ i).type = 'shiny:input';
      }
       Shiny.bindAll();
    }
  });
});
