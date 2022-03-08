$( document ).ready(function() {

  $(document).on('shiny:inputchanged', function(event) {
    if (event.name === 'AdaptativeSelection_ui_1-filesChIP' && event.value) {
      prevText = event.value;
       Shiny.unbindAll();
      //console.log(arr.values(event.value));
      for (let i = 0; i < event.value.length; i++){
        $('#placeholder1').append(
        '<div class=\"form-group shiny-input-container\">' +
          '<label for=\"AdaptativeSelection_ui_1-txt' + i + '\">' + event.value[i].name +  '</label>' +
          '<input id=\"AdaptativeSelection_ui_1-txt' + i + '\" type=\"text\"' +
                  'class=\"form-control\" value=\"\">' + 
        '</div>'
        );
        var pwmname = "";
        $('#AdaptativeSelection_ui_1-txt'+ i).val(pwmname);
        $('#AdaptativeSelection_ui_1-txt'+ i).type = 'shiny:input';
      }
       Shiny.bindAll();
    }
  });
});
