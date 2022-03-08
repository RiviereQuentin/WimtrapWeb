$( document ).ready(function() {

    $(document).on('shiny:value', function(event) {
      
     $('input[name = "checkc"]').on('change', function() {
     $('input[name="checkc"]').not(this).prop('checked', false);
     });
     
    });
});
